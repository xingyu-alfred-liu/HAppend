import numpy as np
from copy import deepcopy


def getSingleMol(supercell, middleSite, bondDict, middleSiteIndex, printMol=True):  # noqa: C901 TODO Fix
    """
    getSingleMol finds a singlemolecule with provided crystal and the, starting
    pymatgen site and index
    Args:
        - printMol: bool
            if True, print out the single molecule information
    Output:
        - dict, key is index of the site, value is the site
    """
    candidates = [middleSite]
    candidatesIndex = [middleSiteIndex]
    tmpSites = []
    singleMol = dict()
    while candidates != []:
        for site in candidates:
            bondCutoffMax = 0
            for pair in bondDict.keys():
                if (bondDict[pair] >= bondCutoffMax) and (str(site.specie) in pair):
                    bondCutoffMax = bondDict[pair]
            allNeighbors = supercell.get_neighbors(site, bondCutoffMax, include_index=True)
            for neighbor in allNeighbors:
                sitePair = (str(site.specie), str(neighbor[0].specie))
                if neighbor[1] <= bondDict[sitePair]:
                    tmpSites += [neighbor]
        # ugly solution
        tmpSiteslist = list()
        for tmpsite in tmpSites:
            if tmpsite not in tmpSiteslist:
                tmpSiteslist += [tmpsite]
        tmpSites = deepcopy(tmpSiteslist)
        for i, site in enumerate(candidates):
            singleMol[candidatesIndex[i]] = site
        candidates = []
        candidatesIndex = []
        for site in tmpSites:
            # check if the site index is in the single molecule index
            if site[2] not in singleMol.keys():
                candidates += [site[0]]
                candidatesIndex += [site[2]]
        if printMol:
            print('Number of atoms found in this molecule:', len(singleMol))
        tmpSites = []
    return singleMol


def getCentralSingleMol(supercell, bondDict, middle=[0.5, 0.5, 0.5], printMol=True):
    """
    Args:
        - supercell: pymatgen Structure
        - bondDict: dict
        - middle: list
        - printMol: bool
            if False, then don't print out the single molecule
            information
    """
    # check if the frac_coord is smaller than zero
    # if smaller than zero, change the middle site coords to negative
    if min(supercell.frac_coords[:, 0]) < 0:
        for i, val in enumerate(middle):
            middle[i] = val * (-1)
    # find site in the middle
    dist = 1
    for i, site in enumerate(supercell.sites):
        if np.linalg.norm(middle-site.frac_coords) < dist:
            dist = np.linalg.norm(middle-site.frac_coords)
            middleSite = site
            middleSiteIndex = i
    if printMol:
        print('The site closest to the middle is', middleSite)
        print('The corresponding site index is', middleSiteIndex)
    # pick up all the atom pairs within bond van der waals distance
    # centralSingleMol is the list of sites which belong to the central molecule
    centralSingleMol = getSingleMol(supercell, middleSite, bondDict, middleSiteIndex, printMol)
    return centralSingleMol


def getBondDict(supercell, bondCutoff):
    bondDict = dict()
    speciesList = list(set(supercell.species))
    for i in range(len(speciesList)):
        speciesList[i] = str(speciesList[i])
    for a in speciesList:
        for b in speciesList:
            if (a, b) in bondCutoff.keys():
                bondDict[(a, b)] = bondCutoff[(a, b)]
    duplicate = dict()
    for pairs in bondDict.keys():
        (a, b) = pairs
        if (b, a) not in bondDict:
            duplicate[(b, a)] = bondDict[pairs]
    bondDict.update(duplicate)
    return bondDict


bondCutoff = {('H', 'H'): 1.2,
              ('H', 'C'): 1.45,
              ('H', 'F'): 1.335,
              ('H', 'N'): 1.375,
              ('H', 'Cl'): 1.475,
              ('H', 'P'): 1.5,
              ('H', 'O'): 1.21,
              ('H', 'S'): 1.5,
              ('H', 'As'): 1.525,
              ('H', 'Br'): 1.525,
              ('H', 'Ga'): 1.535,
              ('H', 'Se'): 1.55,
              ('H', 'I'): 1.59,
              ('H', 'Si'): 1.65,
              ('C', 'C'): 1.7,
              ('C', 'F'): 1.585,
              ('C', 'N'): 1.625,
              ('C', 'O'): 1.53,
              ('C', 'Cl'): 1.79,
              ('C', 'P'): 1.86,
              ('C', 'S'): 2.0,
              ('C', 'As'): 1.775,
              ('C', 'Br'): 1.775,
              ('C', 'Ga'): 1.785,
              ('C', 'Se'): 1.8,
              ('C', 'I'): 1.84,
              ('C', 'Si'): 2.0,
              ('F', 'F'): 1.47,
              ('F', 'N'): 1.51,
              ('F', 'Cl'): 1.61,
              ('F', 'P'): 1.635,
              ('F', 'S'): 1.635,
              ('F', 'As'): 1.66,
              ('F', 'Br'): 1.66,
              ('F', 'Ga'): 1.67,
              ('F', 'Se'): 1.685,
              ('F', 'I'): 1.725,
              ('F', 'Si'): 1.785,
              ('N', 'N'): 1.55,
              ('N', 'Cl'): 1.65,
              ('N', 'P'): 1.675,
              ('N', 'O'): 1.50,
              ('N', 'S'): 1.71,
              ('N', 'As'): 1.7,
              ('N', 'Br'): 1.7,
              ('N', 'Ga'): 1.71,
              ('N', 'Se'): 1.725,
              ('N', 'I'): 1.765,
              ('N', 'Si'): 1.825,
              ('Cl', 'Cl'): 1.75,
              ('Cl', 'P'): 1.775,
              ('Cl', 'S'): 1.775,
              ('Cl', 'As'): 1.8,
              ('Cl', 'Br'): 1.8,
              ('Cl', 'Ga'): 1.81,
              ('Cl', 'Se'): 1.825,
              ('Cl', 'I'): 1.865,
              ('Cl', 'Si'): 1.925,
              ('P', 'P'): 2.25,
              ('P', 'O'): 1.62,
              ('P', 'S'): 2.12,
              ('P', 'As'): 1.825,
              ('P', 'Br'): 1.825,
              ('P', 'Ga'): 1.835,
              ('P', 'Se'): 1.85,
              ('P', 'I'): 1.89,
              ('P', 'Si'): 1.95,
              ('S', 'O'): 1.52,
              ('S', 'S'): 2.11,
              ('S', 'As'): 1.825,
              ('S', 'Br'): 1.825,
              ('S', 'Ga'): 1.835,
              ('S', 'Se'): 1.85,
              ('S', 'I'): 1.89,
              ('S', 'Si'): 1.95,
              ('As', 'As'): 1.85,
              ('As', 'Br'): 1.85,
              ('As', 'Ga'): 1.86,
              ('As', 'Se'): 1.875,
              ('As', 'I'): 1.915,
              ('As', 'Si'): 1.975,
              ('Br', 'Br'): 1.85,
              ('Br', 'Ga'): 1.86,
              ('Br', 'Se'): 1.875,
              ('Br', 'I'): 1.915,
              ('Br', 'Si'): 1.975,
              ('Ga', 'Ga'): 1.87,
              ('Ga', 'Se'): 1.885,
              ('Ga', 'I'): 1.925,
              ('Ga', 'Si'): 1.985,
              ('Se', 'Se'): 1.9,
              ('Se', 'I'): 1.94,
              ('Se', 'Si'): 2.0,
              ('I', 'I'): 1.98,
              ('I', 'Si'): 2.04,
              ('Si', 'Si'): 2.4,
              ('Si', 'O'): 1.69,
              ('O', 'O'): 2.41}
