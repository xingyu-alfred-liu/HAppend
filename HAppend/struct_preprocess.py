import os
import copy
import numpy as np
import math
import json
from ase.io import read
from pymatgen.core import Structure, Molecule, Site
# from pymatgen.io.cif import CifParser
# from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.core.periodic_table import Element, DummySpecie
from pymatgen.analysis.structure_matcher import StructureMatcher, SpeciesComparator
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

from pah101tools.utils import bondCutoff
from pah101tools.utils import getBondDict, getCentralSingleMol, getSingleMol

from rdkit import Chem

def get_singlemolecule(struct):
    superStruct = struct.copy()
    superStruct.make_supercell([3, 3, 3])
    bond_dict = getBondDict(superStruct, bondCutoff)
    single_mol = getCentralSingleMol(superStruct, bond_dict, printMol=False)
    species = []
    position = []
    for i, site in single_mol.items():
        species.append(str(site.specie))
        position.append(site.coords)
    outmol = Molecule(species, position)
    return outmol

def calNormalVector(p1, p2, p3):
    # generate the vector
    vector = [0, 0, 0]
    vector[0] = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
    vector[1] = (p2[2] - p1[2]) * (p3[0] - p1[0]) - (p2[0] - p1[0]) * (p3[2] - p1[2])
    vector[2] = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
    # normalize it
    sigma = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
    sigma = math.sqrt(sigma)
    for i in range(3):
        vector[i] = vector[i] / sigma
    return np.array(vector).reshape(3, )

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


class StructurePreprocess(object):
    """
    Notes:
        1. self.new_struct_dict should keep the format,
           key: original db id
           value: pymatgen structure
        2. adding Hydrogen function is not rubust,
           probably need to make a function myself
    """

    def __init__(self):
        """
        initializing the process
        """

        self.bad_struct_dict = dict()
        self.manual_check_list = list()
        self.add_H_without_origin = list()

        self.db = None
        self.matml_db = None
        self.dbname = None
        self.new_struct_dict = None
        self._duplicate_names = None
        self._radical_system_names = None
        self._invalid_key_doublecheck = None
        self._invalid_mol_dict = None
        self._cocrystal_names = None
        self._still_wrong_after_add_H = None
        self._valid_mol_dict = None
        self._added_H_sites_dict = None

    def load_new_structures(self, path, dbname):
        """
        read in the new structures
        Args:
            path: string
                path to all new structures, should all be
                cif files
            dbname: str
                    the name of the original database, e.g. "CSD", "OMDB"
        TODO: error case handling
        """

        self.dbname = dbname
        self.json_path = path
        self.new_struct_dict = dict()
        for name in os.listdir(path):
            if name.endswith('.json'):
                with open(os.path.join(path, name), "r") as f:
                    data = json.load(f)
                self.new_struct_dict[name.split(".")[0]] = Structure.from_dict(data["geometry"]["wrong_crystal"])

    def load_matml_database(self, db):
        """
        load the current structures and compare one by one
        path is a temporary debugging tool when I work in local

        Args:
            - db: MatMLDatabase
                The database to load structures into (will be wiped clean!)
        """
        self.db = db
        self.matml_db = {}

        # TODO: Currently supports MongoDB only, refactor to support MatMLDatabase
        for _id in self.db.db._collection.find().distinct("_id"):
            self.matml_db[_id] = self.db.db._collection.find_one({"_id": _id})
        print("Loaded %d documents from matml db." % (len(self.matml_db.keys())))

    def move_bad_from_new(self, bad_names=[]):
        """
        a common function to move in general bad structures away from the
        self.new_struct_dict, put into self.bad_struct_dict

        Args:
            bad_names: list
                with the original db struct ids
        """
        if bad_names != []:
            for name in bad_names:
                self.bad_struct_dict[name] = self.new_struct_dict[name]
                self.new_struct_dict.pop(name)

    def chemical_selection(self, chem_group=[
        Element('H'),
        Element('C'),
        Element('N'),
        Element('O'),
        Element('Si'),
        Element('P'),
        Element('S')
    ]
                           ):
        """
        find out structs with elements not of interest

        Args:
            chem_group: list
            A list of elements to be considered as interesting candidates
        """
        self._extra_chem_names = list()
        # sometimes it has D instead of H
        for id, struct in self.new_struct_dict.items():
            if DummySpecie('D') in struct.species:
                struct.replace_species({DummySpecie('D'): Element('H')})
            for ele in list(set(struct.species)):
                if ele not in chem_group:
                    self._extra_chem_names.append(id)
                    break
        # put bad structures into self.bad_struct_dict()
        self.move_bad_from_new(bad_names=self._extra_chem_names)
        print("Found %d structuers contain elements not of interest." % (len(self._extra_chem_names)))
        return self._extra_chem_names

    def get_radical_system(self):
        """
        identify the radical systems

        Args:
            None
        Output:
            name list of radical systems
        """
        self._radical_system_names = list()
        for id, struct in self.new_struct_dict.items():
            super_cell = struct.copy()
            super_cell.make_supercell([5, 5, 5])
            bond_dict = getBondDict(struct, bondCutoff)
            single_mol = getCentralSingleMol(super_cell, bond_dict, printMol=False)
            species = []
            position = []
            delindices = []
            for i, site in single_mol.items():
                species.append(str(site.specie))
                position.append(site.coords)
                delindices.append(i)
            outmol = Molecule(species, position)
            if outmol.nelectrons % 2 != 0:
                self._radical_system_names.append(id)
        # put bad structures into self.bad_struct_dict()
        self.move_bad_from_new(bad_names=self._radical_system_names)
        print("Found %d structuers which are radical systems." % (len(self._radical_system_names)))
        return self._radical_system_names

    def duplicate_check(self):  # noqa: C901 TODO Fix
        """
        duplicate check should also be done with
        existing structures in the matml database
        Remove the hydrogens and then compare?

        1. compare among the new structures
        2. compare the structures one by one
           with the matmldb documents
        """
        self._duplicate_names = []
        _duplicate_names = []
        duplicate_checker = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=True,
            scale=False,
            comparator=SpeciesComparator(),
            ignored_species=['H']
        )
        # first check among new structures
        _struct_id_list = list(self.new_struct_dict.keys())
        for i, s1_id in enumerate(_struct_id_list):
            specie_list1 = list(set(self.new_struct_dict[s1_id].species))
            specie_list1.sort()
            for j in range(i + 1, len(_struct_id_list)):
                s2_id = _struct_id_list[j]
                # check if they have the same species
                specie_list2 = list(set(self.new_struct_dict[s2_id].species))
                specie_list2.sort()
                # if two structure have different species, no need to check
                if not specie_list1 == specie_list2:
                    continue
                if duplicate_checker.fit(self.new_struct_dict[s1_id],
                                         self.new_struct_dict[s2_id]):
                    _duplicate_names.append([s1_id, s2_id])
        for id_pair in _duplicate_names:
            if self.new_struct_dict[id_pair[0]].num_sites >= self.new_struct_dict[id_pair[1]].num_sites:
                self._duplicate_names.append(id_pair[1])
            else:
                self._duplicate_names.append(id_pair[0])
        # now compare with the existing structures
        for s1_id, s1 in self.new_struct_dict.items():
            specie_list1 = list(set(s1.species))
            specie_list1.sort()
            for s2_id, s2_doc in self.matml_db.items():
                s2 = Structure.from_dict(s2_doc["geometry"]["unrelaxed_crystal"])
                specie_list2 = list(set(s2.species))
                specie_list2.sort()
                # check if they have the same species
                if not specie_list1 == specie_list2:
                    continue
                if duplicate_checker.fit(s1, s2):
                    self._duplicate_names.append(s1_id)
        # move the duplicates to bad structures
        self._duplicate_names = list(set(self._duplicate_names))
        self.move_bad_from_new(bad_names=self._duplicate_names)
        print("Found %d duplicated structures." % (len(self._duplicate_names)))
        return self._duplicate_names

    def cocrystal_check(self, fragment_num=8):  # noqa: C901 TODO Fix
        """
        find out cocrystals and put them into bad names
        how: find out the complete unique fragments and
        see if they are the same

        Args:
            - fragment_num: int
                number of fragments closest to the
                center of the supercell
        """
        self._cocrystal_names = []
        molecule_matcher = MoleculeMatcher(tolerance=10)
        for id, struct in self.new_struct_dict.items():
            print("Now checking the ID:", id)
            fragment_list = []
            super_cell = struct.copy()
            super_cell.make_supercell([5, 5, 5])
            bond_dict = getBondDict(struct, bondCutoff)
            # find the eight fragments closest to the center, should
            # already include all possible fragments
            for _ in range(fragment_num):
                single_mol = getCentralSingleMol(super_cell, bond_dict, printMol=False)
                species = []
                position = []
                delindices = []
                for i, site in single_mol.items():
                    species.append(str(site.specie))
                    position.append(site.coords)
                    delindices.append(i)
                outmol = Molecule(species, position)
                fragment_list.append(outmol)
                super_cell.remove_sites(delindices)
            # fragment_list[0].to(fmt='xyz', filename=id+'.xyz')
            # for j, frag in enumerate(fragment_list):
            # frag.to(fmt='xyz', filename=str(j)+'.xyz')
            # after get enough complete fragments, check if there exist single H
            for f in fragment_list:
                if f.num_sites == 1 and str(f.sites[0].specie) == 'H':
                    # if there exist a single H considered as one frag, remove all H
                    # and add the H later
                    self.add_H_without_origin.append(id)
                    self.new_struct_dict[id].remove_species(['H'])
                    break
            if id in self.add_H_without_origin:
                continue
            # 1. check if each has the same numbe of sites
            for i in range(0, len(fragment_list) - 1):
                if fragment_list[i].num_sites != fragment_list[i + 1].num_sites:
                    self._cocrystal_names.append(id)
                    break
            if id in self._cocrystal_names:
                continue
            # 2. if they have the same number of sites, check if the same
            for i, mol1 in enumerate(fragment_list):
                if id in self._cocrystal_names:
                    break
                for j in range(i + 1, len(fragment_list)):
                    mol2 = fragment_list[j]
                    if molecule_matcher.fit(mol1, mol2):
                        continue
                    elif mol1.nelectrons == mol2.nelectrons:
                        print(id, "might be co-crystal, check manually.")
                        self.manual_check_list.append(id)
                        continue
                    else:
                        self._cocrystal_names.append(id)
                        break
            if id not in self._cocrystal_names and struct.num_sites % fragment_list[0].num_sites != 0:
                self.add_H_without_origin.append(id)
                self.new_struct_dict[id].remove_species(['H'])
        # move the cocrystals to bad structures
        print('H are removed in', self.add_H_without_origin)
        self.move_bad_from_new(bad_names=self._cocrystal_names)
        print("Found %d co-crystals." % (len(self._cocrystal_names)))
        return self._cocrystal_names

    def valid_check(self, out=False, double_check=False, invalidFragDict=False):  # noqa: C901 TODO Fix
        """
        check if this material is correct
        use rdkit to see if this can be loaded

        Args:
            - out: bool
                decide if should return the self._invalid_mol_dict
            - doublecheck: bool
                if checking the after-treated mol instead of
                pulling one mol from the struct
        """
        if invalidFragDict is False:
            if double_check:
                check_mol_dict = copy.deepcopy(self._valid_mol_dict)
                self._invalid_key_doublecheck = []
            else:
                self._invalid_mol_dict = dict()
                check_mol_dict = dict()
                for key, struct in self.new_struct_dict.items():
                    pmgmol = get_singlemolecule(struct)
                    check_mol_dict[key] = pmgmol
        else:
            _invalid_mol_dict = dict()
            check_mol_dict = copy.deepcopy(invalidFragDict)
        for key, pmgmol in check_mol_dict.items():
            # find the nitro group
            nitroIDs = []
            # pmgmol.to(fmt='xyz', filename='tmptmp.xyz')
            if 'N' in pmgmol.symbol_set and 'O' in pmgmol.symbol_set:
                # decide if it is nitro group
                Nidlist = []
                for i, site in enumerate(pmgmol.sites):
                    if site.specie == Element("N"):
                        Nidlist.append(i)
                if not Element('H') in pmgmol.species:
                    tmpmol = pmgmol.copy()
                    tmpmol.append('H', [100, 0.1, 200])
                    bond_dict = getBondDict(tmpmol, bondCutoff)
                else:
                    bond_dict = getBondDict(pmgmol, bondCutoff)
                max_key = max(bond_dict.keys(), key=(lambda k: bond_dict[k]))
                bond_length_max = bond_dict[max_key]
                for Nid in Nidlist:
                    N_neighbors = self.get_bonded_neighbors(pmgmol, bond_dict, Nid, bond_length_max)
                    if len(N_neighbors) != 3:
                        continue
                    bonded_O = []
                    for neighbor in N_neighbors:
                        if neighbor.specie == Element("O"):
                            bonded_O.append(neighbor)
                    if len(bonded_O) != 2:
                        continue
                    is_nitro = True
                    tmpOids = []
                    for each_O in bonded_O:
                        for j, site in enumerate(pmgmol.sites):
                            if each_O == site:
                                Oid = j
                                tmpOids.append(Oid)
                                break
                        Oneighbors = self.get_bonded_neighbors(pmgmol, bond_dict, Oid, bond_length_max)
                        if len(Oneighbors) != 1 or Oneighbors[0].specie != Element("N"):
                            is_nitro = False
                            break
                    if is_nitro:
                        tmpnitroIDs = [Nid]
                        tmpnitroIDs += tmpOids
                        if len(tmpnitroIDs) == 3:
                            nitroIDs += tmpnitroIDs
            obadaptor = BabelMolAdaptor(pmgmol)
            # TODO: find out a method interfacing between openbabel and rdkit
            obadaptor.write_file('tmp.mol', file_format='mol')
            rdkmol = Chem.MolFromMolFile('tmp.mol', removeHs=False)
            os.remove('tmp.mol')
            if rdkmol is None:
                # rdkit cannot read this structure
                self.manual_check_list.append(key)
                try:
                    self.bad_struct_dict[key] = self.new_struct_dict[key]
                    del self.new_struct_dict[key]
                    continue
                except:  # noqa: E722 TODO Fix
                    continue
            wrongAtomIdx = []
            wrongAtomHybrid = []
            wrongAtomValence = []
            for i in range(rdkmol.GetNumAtoms()):
                if i in nitroIDs:
                    continue
                rdkatom = rdkmol.GetAtomWithIdx(i)
                if rdkatom.GetSymbol() == "C":
                    if rdkatom.GetExplicitValence() == 4:
                        pass
                    else:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
                if rdkatom.GetSymbol() == "N":
                    if rdkatom.GetExplicitValence() != 3 and rdkatom.GetExplicitValence() != 5:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
                if rdkatom.GetSymbol() == "O":
                    if rdkatom.GetExplicitValence() == 2:
                        pass
                    else:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
                if rdkatom.GetSymbol() == "Si":
                    if rdkatom.GetExplicitValence() == 4:
                        pass
                    else:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
                if rdkatom.GetSymbol() == "P":
                    if rdkatom.GetExplicitValence() == 3 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        pass
                    elif rdkatom.GetExplicitValence() == 4 and \
                            rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3D:
                        pass
                    elif rdkatom.GetExplicitValence() == 5 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        pass
                    else:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
                if rdkatom.GetSymbol() == "S":
                    if rdkatom.GetExplicitValence() == 2 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                        pass
                    elif rdkatom.GetExplicitValence() == 2 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        pass
                    elif rdkatom.GetExplicitValence() == 4 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        pass
                    elif rdkatom.GetExplicitValence() == 6 and rdkatom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        pass
                    else:
                        wrongAtomIdx.append(i)
                        wrongAtomHybrid.append(rdkatom.GetHybridization())
                        wrongAtomValence.append(rdkatom.GetExplicitValence())
            if len(wrongAtomIdx) != 0:
                if invalidFragDict is False:
                    if not double_check:
                        self._invalid_mol_dict[key] = dict()
                        self._invalid_mol_dict[key]["pmgmol"] = pmgmol
                        self._invalid_mol_dict[key]["wrongID"] = wrongAtomIdx
                        self._invalid_mol_dict[key]["hybrid"] = wrongAtomHybrid
                        self._invalid_mol_dict[key]["valence"] = wrongAtomValence
                    else:
                        self._invalid_key_doublecheck.append(key)
                else:
                    _invalid_mol_dict[key] = dict()
                    _invalid_mol_dict[key]["pmgmol"] = pmgmol
                    _invalid_mol_dict[key]["wrongID"] = wrongAtomIdx
                    _invalid_mol_dict[key]["hybrid"] = wrongAtomHybrid
                    _invalid_mol_dict[key]["valence"] = wrongAtomValence
        if double_check:
            print("After double check, found %d structs are still wrong:"
                  % (len(self._invalid_key_doublecheck)))
            print(self._invalid_key_doublecheck)
        else:
            if invalidFragDict is False:
                print("Found %d possibly wrong molecules."
                      % (len(self._invalid_mol_dict.keys())))
        if double_check:
            # move the wrong structs to the bad list, also delete the key from
            # self._added_H_sites_dict
            try:
                self.move_bad_from_new(bad_names=self._invalid_key_doublecheck)
            except:  # noqa: E722 TODO Fix
                pass
            for structID in self._invalid_key_doublecheck:
                del self._added_H_sites_dict[structID]
            if out:
                return self._invalid_key_doublecheck
        else:
            if out:
                if invalidFragDict is False:
                    return self._invalid_mol_dict
                else:
                    return _invalid_mol_dict

    def get_all_fragments(self, struct, supercellsize):
        """
        Args:
            struct: pymatgen Structure
                looking for complete fragment from this struct

        Return:
            list of fragments
        """
        fragmentListAll = []
        struct_super = struct.copy()
        translate = np.sum(struct.lattice.matrix, 0)
        struct_super.make_supercell(supercellsize)
        bond_dict = getBondDict(struct_super, bondCutoff)
        frag = get_singlemolecule(struct)
        if supercellsize == [3, 3, 3]:
            translate = np.sum(struct.lattice.matrix, 0)
        elif supercellsize == [5, 5, 5]:
            translate = np.sum(struct.lattice.matrix, 0) * 2
        elif supercellsize == [7, 7, 7]:
            translate = np.sum(struct.lattice.matrix, 0) * 3
        elif supercellsize == [9, 9, 9]:
            translate = np.sum(struct.lattice.matrix, 0) * 4
        while struct_super.num_sites != 0:
            middleSite = struct_super.sites[0]
            middleSiteIndex = 0
            new_frag = getSingleMol(struct_super, middleSite, bond_dict, middleSiteIndex, printMol=False)
            species = []
            position = []
            for k, site in new_frag.items():
                species.append(str(site.specie))
                position.append(site.coords - translate)
            struct_super.remove_sites(new_frag.keys())
            outmol = Molecule(species, position)
            if outmol.num_sites == frag.num_sites:
                fragmentListAll.append(outmol)
            else:
                pass
        print("Get %d complete fragments." % (len(fragmentListAll)))
        return fragmentListAll

        # fragmentListAll = []
        # struct_super = struct.copy()
        # translate = np.sum(struct.lattice.matrix, 0)
        # struct_super.make_supercell([3, 3, 3])
        # bond_dict = getBondDict(struct_super, bondCutoff)
        # for i, originSite in enumerate(struct.sites):
        #     translatedSiteCoords = originSite.coords + translate
        #     for j, superSite in enumerate(struct_super.sites):
        #         if np.linalg.norm(translatedSiteCoords-superSite.coords) < 0.1:
        #             middleSite = superSite
        #             middleSiteIndex = j
        #             break
        #     single_mol = getSingleMol(struct_super, middleSite, bond_dict, \
        #                     middleSiteIndex, printMol=False)
        #     species = []
        #     position = []
        #     for k, site in single_mol.items():
        #         species.append(str(site.specie))
        #         position.append(site.coords-translate)
        #     # struct_super.remove_sites(single_mol.keys())
        #     outmol = Molecule(species, position)
        #     fragmentListAll.append(outmol)
        # print("Found %d complete fragments." %(len(fragmentListAll)))
        # fragList = []
        # for mol in fragmentListAll:
        #     exist=False
        #     for frag in fragList:
        #         if np.linalg.norm(frag.center_of_mass-mol.center_of_mass) < 0.01:
        #             exist = True
        #             break
        #     if not exist:
        #         fragList.append(mol)
        # tmpmolList = []
        # for x in [-1, 0, 1]:
        #     for y in [-1, 0, 1]:
        #         for z in [-1, 0, 1]:
        #             if [x,y,z] == [0,0,0]:
        #                 continue
        #             transVec = np.sum(struct.lattice.matrix * [x, y, z], 0)
        #             transVec = np.concatenate((transVec, np.array([1]).reshape(1,)), 0)
        #             rotationMatrix = np.zeros((4, 3), int)
        #             np.fill_diagonal(rotationMatrix, 1)
        #             affineMatrix = np.concatenate((rotationMatrix, transVec.reshape(-1, 1)), 1)
        #             symOp = SymmOp(affineMatrix)
        #             for frag in fragList:
        #                 species = []
        #                 coords = []
        #                 for site in frag.sites:
        #                     species.append(site.specie)
        #                     coords.append(symOp.operate(site.coords))
        #                 tmpmolList.append(Molecule(species, coords))
        # for fragx in tmpmolList[:]:
        #     for fragy in fragList:
        #         if np.linalg.norm(fragx.center_of_mass-fragy.center_of_mass) < 1:
        #             tmpmolList.remove(fragx)
        # dellist = []
        # for i, fragx in enumerate(tmpmolList):
        #     for j in range(i+1, len(tmpmolList)):
        #         fragy = tmpmolList[j]
        #         if np.linalg.norm(fragx.center_of_mass-fragy.center_of_mass) < 1:
        #             dellist.append(fragy)
        # for frag in dellist:
        #     try:
        #         tmpmolList.remove(frag)
        #     except:
        #         pass
        # fragList += tmpmolList
        # return fragList

    def get_bonded_neighbors(self, pmgmol, bond_dict, siteID, bond_length_max):
        bonded_neighbor = []
        wrongSite = pmgmol.sites[siteID]
        neighbors = pmgmol.get_neighbors(wrongSite, bond_length_max)
        # delete the wrong neighbors if they are not bonded with wrongSite
        for j, neighbor in enumerate(neighbors):
            sitePair = (str(wrongSite.specie), str(neighbor[0].specie))
            if neighbor[1] <= bond_dict[sitePair]:
                bonded_neighbor.append(neighbor[0])
        return bonded_neighbor

    @staticmethod  # noqa: C901 TODO Fix
    def add_Hydrogen_Site(wrongSite, hybrid, neighbors, hydrogenBondLength, key, siteID, com, mol, bond_dict):
        """
        Args:
            wrongSite: pymatgen Site
                the site which needs hydrogen
            hybrid: rdkit.Chem.HybridizationType
                hybridization type of the wrong site
            neighbors: list
                list of neighbors output from
                pymatgen.structure.get_neighbors
            hydrogenBondLength: flot
                value of the wrongSite-H bond threshold
                in dbaAutomator
            key: str
                the ID of the material
            siteID:str
                ID of the wrongSite
            com: np.array
                center of mass
            mol: pymatgen.Molecule
                wrong molecule
            bond_dict: dict
                bond length threshold of each bonds

        TODO:
            1. implement conditions where P has 5 explicit valences
        """
        hydrogen_site_list = []
        if wrongSite.specie in [Element("C"), Element("Si")]:
            if hybrid == Chem.rdchem.HybridizationType.SP:
                if len(neighbors) != 1:
                    print("Check %s for its %d site." % (key, siteID))
                else:
                    dist = np.linalg.norm((wrongSite.coords - neighbors[0][0].coords), 2)
                    coord = wrongSite.coords - 0.7 * hydrogenBondLength / dist * (
                            neighbors[0][0].coords - wrongSite.coords)
                    hydrogen_site_list.append(Site("H", coord))
            elif hybrid == Chem.rdchem.HybridizationType.SP2:
                if len(neighbors) == 0:
                    print("Check %s for its %d site." % (key, siteID))
                elif len(neighbors) == 1:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerCcoord = wrongSite.coords - 0.3 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    normVec = np.cross(distVec, com)
                    normVec = normVec / np.linalg.norm(normVec)
                    Hcoord_1 = centerCcoord + 0.4 * normVec * hydrogenBondLength
                    # locate the fist H atom
                    # Hcoord_1 = wrongSite.coords + normVec*0.8*hydrogenBondLength
                    # find out the next two H atoms by rotating the first one
                    Hcoord_2 = np.dot(rotation_matrix(distVec, math.pi), Hcoord_1 - wrongSite.coords) + wrongSite.coords
                    hydrogen_site_list.append(Site("H", Hcoord_1))
                    hydrogen_site_list.append(Site("H", Hcoord_2))
                elif len(neighbors) == 2:
                    middleCoord = 0.5 * (neighbors[0][0].coords + neighbors[1][0].coords)
                    dist = np.linalg.norm((wrongSite.coords - middleCoord), 2)
                    coord = wrongSite.coords - 0.6 * hydrogenBondLength / dist * (middleCoord - wrongSite.coords)
                    hydrogen_site_list.append(Site("H", coord))
            elif hybrid == Chem.rdchem.HybridizationType.SP3:
                if len(neighbors) == 0:
                    print("Check %s for its %d site." % (key, siteID))
                elif len(neighbors) == 1:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerNcoord = wrongSite.coords - 0.3 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    normVec = np.cross(distVec, wrongSite.coords - com)
                    normVec = normVec / np.linalg.norm(normVec)
                    Hcoord_1 = centerNcoord + 0.55 * normVec * hydrogenBondLength
                    # find out the next two H atoms by rotating the first one
                    Hcoord_2 = np.dot(rotation_matrix(distVec, math.pi * 2 / 3),
                                      Hcoord_1 - wrongSite.coords) + wrongSite.coords
                    Hcoord_3 = np.dot(rotation_matrix(distVec, math.pi * 4 / 3),
                                      Hcoord_1 - wrongSite.coords) + wrongSite.coords
                    # add the to the hydrogen list
                    hydrogen_site_list.append(Site("H", Hcoord_1))
                    hydrogen_site_list.append(Site("H", Hcoord_2))
                    hydrogen_site_list.append(Site("H", Hcoord_3))
                elif len(neighbors) == 2:
                    middleCoord = 0.5 * (neighbors[0][0].coords + neighbors[1][0].coords)
                    dist = np.linalg.norm((wrongSite.coords - middleCoord), 2)
                    HmiddleCoord = wrongSite.coords - 0.4 * hydrogenBondLength / dist * (middleCoord - wrongSite.coords)
                    normalVec = calNormalVector(wrongSite.coords, neighbors[0][0].coords, neighbors[1][0].coords)
                    Hcoord_1 = HmiddleCoord + normalVec * 0.5 * hydrogenBondLength
                    Hcoord_2 = HmiddleCoord - normalVec * 0.5 * hydrogenBondLength
                    hydrogen_site_list.append(Site("H", Hcoord_1))
                    hydrogen_site_list.append(Site("H", Hcoord_2))
                elif len(neighbors) == 3:
                    middleCoord = (neighbors[0][0].coords + neighbors[1][0].coords + neighbors[2][0].coords) / 3
                    dist = np.linalg.norm((wrongSite.coords - middleCoord), 2)
                    coord = wrongSite.coords - 0.6 * hydrogenBondLength / dist * (middleCoord - wrongSite.coords)
                    hydrogen_site_list.append(Site("H", coord))
        if wrongSite.specie in [Element("O"), Element("S")]:
            if hybrid == Chem.rdchem.HybridizationType.SP2:
                # find out the existing bond and rotate 110 degrees to get the H coord
                if len(neighbors) != 1:
                    print("Check %s for its %d site." % (key, siteID))
                else:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerOcoord = wrongSite.coords - 0.2 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    # normVec = np.cross(distVec, wrongSite.coords-com)
                    # normVec = normVec/np.linalg.norm(normVec)
                    # coord = centerOcoord + 0.6*normVec*hydrogenBondLength
                    normVec = wrongSite.coords - com
                    normVec = normVec / np.linalg.norm(normVec)
                    coord = centerOcoord + 0.6 * normVec * hydrogenBondLength
                    hydrogen_site_list.append(Site("H", coord))
            elif hybrid == Chem.rdchem.HybridizationType.SP3:
                # find out the existing bond and rotate 110 degrees to get the H coord
                if len(neighbors) != 1:
                    print("Check %s for its %d site." % (key, siteID))
                else:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerOcoord = wrongSite.coords - 0.2 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    # normVec = np.cross(distVec, wrongSite.coords-com)
                    # normVec = normVec/np.linalg.norm(normVec)
                    # coord = centerOcoord + 0.6*normVec*hydrogenBondLength
                    normVec = wrongSite.coords - com
                    normVec = normVec / np.linalg.norm(normVec)
                    coord = centerOcoord + 0.6 * normVec * hydrogenBondLength
                    hydrogen_site_list.append(Site("H", coord))
            else:
                print("Check %s for its %d site." % (key, siteID))
        if wrongSite.specie in [Element("N"), Element("P")]:
            if hybrid == Chem.rdchem.HybridizationType.SP3:
                if len(neighbors) == 0:
                    print("Check %s for its %d site." % (key, siteID))
                elif len(neighbors) == 1:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerNcoord = wrongSite.coords - 0.3 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    normVec = np.cross(distVec, com)
                    normVec = normVec / np.linalg.norm(normVec)
                    # locate the fist H atom
                    Hcoord_1 = centerNcoord + 0.5 * normVec * hydrogenBondLength
                    Hcoord_2 = centerNcoord - 0.5 * normVec * hydrogenBondLength
                    hydrogen_site_list.append(Site("H", Hcoord_1))
                    hydrogen_site_list.append(Site("H", Hcoord_2))
                elif len(neighbors) == 2:
                    middleNoord = 0.5 * (neighbors[0][0].coords + neighbors[1][0].coords)
                    dist = np.linalg.norm((wrongSite.coords - middleNoord), 2)
                    coord = wrongSite.coords - 0.6 * hydrogenBondLength / dist * (middleNoord - wrongSite.coords)
                    hydrogen_site_list.append(Site("H", coord))
            elif hybrid == Chem.rdchem.HybridizationType.SP2:
                if len(neighbors) == 0:
                    print("Check %s for its %d site." % (key, siteID))
                elif len(neighbors) == 1:
                    # decide another vector to build the normal vector
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    centerNcoord = wrongSite.coords - 0.3 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    # find out a normal vector perpendicular to worngsite-neighbor
                    max_key = max(bond_dict.keys(), key=(lambda k: bond_dict[k]))
                    bond_length_max = bond_dict[max_key]
                    tmpNeighbors = mol.get_neighbors(neighbors[0][0], bond_length_max)
                    bonded_neighbor = []
                    for _, tmpneighbor in enumerate(tmpNeighbors):
                        sitePair = (str(neighbors[0][0].specie), str(tmpneighbor[0].specie))
                        if tmpneighbor[1] <= bond_dict[sitePair]:
                            bonded_neighbor.append(tmpneighbor[0])
                    # delete wrongSite from it
                    for q, tmpsite in enumerate(bonded_neighbor):
                        if np.linalg.norm(tmpsite.coords - wrongSite.coords) < 0.01:
                            del bonded_neighbor[q]
                            break
                    if len(bonded_neighbor) >= 2:
                        tmpVec = calNormalVector(neighbors[0][0].coords, bonded_neighbor[0].coords,
                                                 bonded_neighbor[1].coords)
                        normVec = np.cross(distVec, tmpVec)
                    else:
                        normVec = np.cross(distVec, wrongSite.coords - com)
                    normVec = normVec / np.linalg.norm(normVec)
                    Hcoord_1 = centerNcoord + 0.4 * normVec * hydrogenBondLength
                    Hcoord_2 = centerNcoord - 0.4 * normVec * hydrogenBondLength
                    hydrogen_site_list.append(Site("H", Hcoord_1))
                    hydrogen_site_list.append(Site("H", Hcoord_2))
                elif len(neighbors) == 2:
                    middleCoord = 0.5 * (neighbors[0][0].coords + neighbors[1][0].coords)
                    dist = np.linalg.norm((wrongSite.coords - middleCoord), 2)
                    coord = wrongSite.coords - 0.7 * hydrogenBondLength / dist * (middleCoord - wrongSite.coords)
                    hydrogen_site_list.append(Site("H", coord))
            elif hybrid == Chem.rdchem.HybridizationType.SP:
                # find out the existing bond and rotate 110 degrees to get the H coord
                if len(neighbors) != 1:
                    print("Check %s for its %d site." % (key, siteID))
                else:
                    distVec = neighbors[0][0].coords - wrongSite.coords
                    coord = wrongSite.coords - 0.7 * hydrogenBondLength / np.linalg.norm(distVec) * distVec
                    hydrogen_site_list.append(Site("H", coord))
        return hydrogen_site_list

    def add_Hydrogen_Molecule(self, invalid_mol_dict=None, out=False, fragDict=False):  # noqa: C901 TODO Fix
        """
        Args:
            invalid_mol_dict: dict
                same data structure as the one of
                self._invalid_mol_dict.
                if is None, use self._invalid_mol_dict
            out: bool
                if return the valid molecule dictionary
        """
        if fragDict is False:
            self._valid_mol_dict = dict()
            self._added_H_sites_dict = dict()
        else:
            _fragDictH = dict()
        if invalid_mol_dict is None:
            invalid_mol_dict = self._invalid_mol_dict
        for key, content in invalid_mol_dict.items():
            pmgmol = content["pmgmol"]
            # pmgmol.to(fmt='xyz', filename='wrong.xyz')
            if not Element('H') in pmgmol.species:
                tmpmol = pmgmol.copy()
                tmpmol.append('H', [100, 0.1, 200])
                bond_dict = getBondDict(tmpmol, bondCutoff)
            else:
                bond_dict = getBondDict(pmgmol, bondCutoff)
            max_key = max(bond_dict.keys(), key=(lambda k: bond_dict[k]))
            bond_length_max = bond_dict[max_key]
            neighbors_dict = dict()
            for siteID in content["wrongID"]:
                bonded_neighbor = self.get_bonded_neighbors(pmgmol, bond_dict, siteID, bond_length_max)
                neighbors_dict[siteID] = bonded_neighbor
            # adding the hydrogens
            hydrogen_site_list = []
            for k, siteID in enumerate(content["wrongID"]):
                # discuss when have different types of hybridization
                wrongSite = pmgmol.sites[siteID]
                sitePair = ('H', str(wrongSite.specie))
                hydrogenBondLength = bond_dict[sitePair]
                # call the add H Site function
                hydrogen_site_list += self.add_Hydrogen_Site(
                    wrongSite=wrongSite,
                    hybrid=content["hybrid"][k],
                    neighbors=neighbors_dict[siteID],
                    hydrogenBondLength=hydrogenBondLength,
                    key=key,
                    siteID=siteID,
                    com=pmgmol.center_of_mass,
                    mol=pmgmol,
                    bond_dict=bond_dict)
            keepstruct = True
            for site in hydrogen_site_list:
                try:
                    pmgmol.append(site.specie, site.coords)
                except:  # noqa: E722 TODO Fix
                    print("Cannot add Hydrogen for:", key)
                    self.manual_check_list.append(key)
                    self.bad_struct_dict[key] = self.new_struct_dict[key]
                    del self.new_struct_dict[key]
                    keepstruct = False
                    break
            # pmgmol.to(fmt='xyz', filename=str(key)+'_new.xyz')
            if fragDict is False:
                if keepstruct:
                    self._added_H_sites_dict[key] = hydrogen_site_list
                    self._valid_mol_dict[key] = pmgmol
            else:
                _fragDictH[key] = hydrogen_site_list
        if fragDict is False:
            if out:
                return self._valid_mol_dict
        else:
            if out:
                return _fragDictH

    def add_Hydrogen_Structure(self, added_H_sites_dict=None):  # noqa: C901 TODO Fix
        """
        Args:
            added_H_sites_dict: dict
                key is the structID, value is list of added H for the central mol
        TODO:
            sometimes [3,3,3] is not enough for all complete molecules, need to use [5,5,5]
            should implement a separate condition
        """
        wrong_list = []
        self._still_wrong_after_add_H = []
        if added_H_sites_dict is None:
            added_H_sites_dict = self._added_H_sites_dict
        # first find out all fragments of the struct
        self.all_fragment_dict = dict()
        for structID in self._added_H_sites_dict.keys():
            center_frag_H_sites = copy.deepcopy(added_H_sites_dict[structID])
            struct = self.new_struct_dict[structID]
            center_frag = get_singlemolecule(struct)
            translate = np.sum(struct.lattice.matrix, 0)
            for i, siteH in enumerate(center_frag_H_sites):
                center_frag_H_sites[i] = Site("H", siteH.coords - translate)
            correct = False
            for cellsize in [3, 5, 7]:  # [3,5,7,9]
                fragmentList = self.get_all_fragments(struct=struct, supercellsize=[cellsize] * 3)
                invalidFragDict = dict()
                for i, frag in enumerate(fragmentList):
                    invalidFragDict[i] = frag
                invalidFragDict = self.valid_check(out=True,
                                                   double_check=False,
                                                   invalidFragDict=invalidFragDict)
                addHdict = self.add_Hydrogen_Molecule(invalid_mol_dict=invalidFragDict,
                                                      out=True,
                                                      fragDict=True)
                tmpStruct = struct.copy()
                appendStruct = struct.copy()
                for _, Hlist in addHdict.items():
                    for eachH in Hlist:
                        tmpStruct.append("H", eachH.coords, coords_are_cartesian=True)
                        add_this_H = True
                        for val in tmpStruct.sites[-1].frac_coords:
                            if val < 0 or val > 1:
                                add_this_H = False
                                break
                        if add_this_H:
                            appendStruct.append("H", eachH.coords, coords_are_cartesian=True)
                print("now have the cif for:", structID)
                # check if this crytal is correct
                appendStruct.to(fmt='cif', filename=self.json_path+"/"+structID + '_' + str(cellsize) + '.cif')
                if appendStruct.num_sites % (len(center_frag_H_sites) + center_frag.num_sites) == 0:
                    self.new_struct_dict[structID] = appendStruct
                    correct = True
                    break
            if not correct:
                wrong_list.append(structID)
                del self.new_struct_dict[structID]
        if len(wrong_list) != 0:
            print("After added Hydrogens, these structures might still be wrong:")
            print(wrong_list)
            self._still_wrong_after_add_H = wrong_list

    def generate_document(self, struct_dict=None):  # noqa: C901 TODO Fix
        """
        Notice, should insert all necessary keys and leave
        values as None, e.g., relaxed geometry should all be None

        For the results part, each value should have "val" and
        "log", where "val" is the key for calculated values,
        "log" is the key for the whole record of the calculation,
        for FHI-aims results is the JSON format results, for BerkeleyGW
        is the whole output file. The "log" is just for debug purpose.
        Use "bandgap" as an example.

        For the predictions, each predicted result should contain the
        prediction and the corresponding model path. Use "bandgap" as an
        example.

        Notice that for the "problems", if the structure is found to be 100%
        correct, this part will not exist

        need to insert:
            - "struct_id": "ABECAL",
            - "origin_db": "CSD",
            - "_id": "sdf34ff24tgdg2s24f" (should be handled directly by MongoDB)
            - "geometry": {
                "unrelaxed_crystal": {
                    pymatgen.Structure.to_dict()
                },
                "unrelaxed_molecule": {
                    pymatgen.Molecule.to_dict()
                },
                "relaxed_crystal": None,
                "relaxed_molecule": None
            },
            - "dft": {
                "bandgap": {
                    "val": 1.5,
                    "log": ["..."]
                },
                "Et": None,
                "DF": None,
                "VBdisp": None,
                "CBdisp": None,
                "hab": None,
                "gap_s": None,
                "Et_s": None,
                "DF_s": None,
                "IP_s": None,
                "EA_s": None,
                "polarisation": None,
                "apc": None,
                "density": None,
                "epsilon": None,
                "weight_s": None,
            },
            - "gwbse": {
                "Es": None,
                "Et": None,
                "DF": None
            },
            - "predictions": {
                "bandgap": {
                    "val": "SMALL",
                    "model": "/PATH/TO/MODEL"
                },
                "gwbse_DF": None
            },
            - "problems": ['missingH', 'cocrystal', ...]
        }
        """
        if struct_dict is None:
            struct_dict = self.new_struct_dict
        # combine all bad names to decide if need to add "problem"
        bad_name_all = \
            self._extra_chem_names + \
            self._radical_system_names + \
            self._cocrystal_names + \
            self._duplicate_names + \
            list(self._invalid_mol_dict.keys()) + \
            self._still_wrong_after_add_H
        bad_name_all = list(set(bad_name_all))
        # do this fomatting for each material
        for key, pmg_struct in struct_dict.items():
            """
                format a standard document
            """
            struct_json = dict()
            # the original db ID
            struct_json["struct_id"] = key
            # dbname, the name of the original db
            struct_json["origin_db"] = self.dbname

            # the 'geometry' dictionary should include unrelaxed
            # crystal and single molecule, together with relaxed
            # crystal and single molecule, at first there should be
            # NO relaxed structures
            struct_json["geometry"] = dict()

            # under "geometry" should be four entries, two unrelaxed
            # crystals, two relaxed single molecules
            # unrelaxed structures
            struct_json["geometry"]["unrelaxed_crystal"] = pmg_struct.as_dict()
            singlemol = get_singlemolecule(pmg_struct)
            struct_json["geometry"]["unrelaxed_molecule"] = singlemol.as_dict()

            # relaxed structures, unavailable for now
            struct_json["geometry"]["relaxed_crystal"] = None
            struct_json["geometry"]["relaxed_molecule"] = None

            # calculated dft properties
            dft_property_list = ["bandgap", "Et", "DF", "VBdisp", "CBdisp", "hab", "gap_s", "Et_s", "DF_s",
                                 "IP_s", "EA_s", "polarisation", "apc", "density", "epsilon", "weight_s"]
            struct_json["dft"] = dict()
            for property in dft_property_list:
                struct_json["dft"][property] = dict()
                struct_json["dft"][property]["val"] = None
                struct_json["dft"][property]["log"] = None

            # calculated gwbse properties
            gwbse_property_list = ["Es", "Et", "DF"]
            struct_json["gwbse"] = dict()
            for property in gwbse_property_list:
                struct_json["gwbse"][property] = dict()
                struct_json["gwbse"][property]["val"] = None
                struct_json["gwbse"][property]["log"] = None
                struct_json["gwbse"][property]["path"] = None

            # predictions
            prediction_list = ["bandgap", "gwbse_DF"]
            struct_json["predictions"] = dict()
            for prediction in prediction_list:
                struct_json["predictions"][prediction] = dict()
                struct_json["predictions"][prediction]["val"] = None
                struct_json["predictions"][prediction]["model"] = None

            # decide if need a "problem" section
            if key in bad_name_all:
                struct_json["problems"] = list()
                if key in self._extra_chem_names:
                    struct_json["problems"].append("extra_chem")
                if key in self._radical_system_names:
                    struct_json["problems"].append("radical_system")
                if key in self._cocrystal_names:
                    struct_json["problems"].append("cocrystal")
                if key in self._duplicate_names:
                    struct_json["problems"].append("duplicate")
                if key in list(self._invalid_mol_dict.keys()):
                    struct_json["problems"].append("invalid_bonds")
                if key in self._still_wrong_after_add_H:
                    struct_json["problems"].append("wrong_after_add_H")

            # put this struct_json under this key
            struct_dict[key] = struct_json

    def save_struct(self, doc=None):
        """
        save the new structures

        Args:
            doc: dict
                if not specified, the self.new_struct_dict will
                be stored
        """
        # if doc is None:
        #     doc = self.new_struct_dict
        # print("Saving these structs...")
        # for name, data in doc.items():
        #     print(name)
        #     self.matml_db._collection.insert_one(data)

        # -- TEMPORARY STORE AS JSON -- #
        # -- delete after mongo is available -- #
        if doc is None:
            doc = self.new_struct_dict
        print("Saving structs...")
        for name, data in doc.items():
            print("saving " + name)
            inserted_acknowledge, inserted_id = self.db.add_structs(data)
            print("Inserted Acknowledgement:", inserted_acknowledge)
            print("Inserted ID:", inserted_id)
