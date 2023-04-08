from pah101tools.struct_preprocess import StructurePreprocess
from pymatgen.core import Element

def call_struct_preprocess(new_struct_path, dbname):
    struct_processor = StructurePreprocess()
    struct_processor.load_new_structures(path=new_struct_path, dbname=dbname)
    
    # 1. chemical configuration screening, only keep
    # C, H others put into bad ones
    chem_group = [Element('H'), Element('C')]
    _ = struct_processor.chemical_selection(chem_group=chem_group)


    # 2. find out the invalid mols
    invalid_mol_dict = struct_processor.valid_check(out=True)
    print("invalid_mol_dict:", invalid_mol_dict.keys())
    print()

    # 3. add the hydrogens for the mols not struct yet
    valid_mol_dict = struct_processor.add_Hydrogen_Molecule(out=True)
    print("valid_mol_dict:", valid_mol_dict.keys())
    print()

    # 4. check if this is a valid mol, valid check doublecheck=True
    struct_processor.valid_check(double_check=True)

    # 5. add the Hydrogens for the crystals
    struct_processor.add_Hydrogen_Structure()

def main():
    # define the path to the jsons
    new_struct_path = "example_json"
    # define the database name
    dbname = "CSD"
    call_struct_preprocess(new_struct_path, dbname)

if __name__ == "__main__":
    main()