# PAH101Plot
Code to reproduce the plots in paper. In addition, the add Hydrogen program is attahced. 

### Installation
Go in the project folder and type:  
```
python setup.py develop
```

### Example
__Plot bandstructure and absorption__  
In the `example` folder, to gather the bandstructure and absorption for ABECAL, type in:  

```
python plot_bandstructure_absorption.py -p ../data -id ABECAL -b True -a True -bp band.png -ap absorp.png
```
Details for each argument can be gathered by:  

```
python plot_bandstructure_absorption.py -h  

usage: plot_bandstructure_absorption.py [-h] [-p ROOTPATH] [-id STRUCT_ID] [-b BANDSTRUCTURE]
                                        [-a ABSORPTION] [-bp BANDSTRUCTURE_PATH] [-ap ABSORPTION_PATH]

optional arguments:
  -h, --help            show this help message and exit
  -p ROOTPATH, --rootpath ROOTPATH
                        root path point to all JSON
  -id STRUCT_ID, --struct_id STRUCT_ID
                        struct_id match with materials name
  -b BANDSTRUCTURE, --bandstructure BANDSTRUCTURE
                        if plot bandstructure
  -a ABSORPTION, --absorption ABSORPTION
                        if plot absorption spectrum
  -bp BANDSTRUCTURE_PATH, --bandstructure_path BANDSTRUCTURE_PATH
                        target path to save bandstructure plot
  -ap ABSORPTION_PATH, --absorption_path ABSORPTION_PATH
                        target path to save absorption spectrum
```
__Add Hydrogen to crystal__  
In the `example` folder, to add the missing Hydrogen to a structure, using the example ABECAL, where all Hydrogen atoms were removed for demo purpose. Type in:

```
python perform_struct_preprocess.py
```
which will load in the `ABECAL_example.json` file, get the wrong structure, try to add the Hydrogen back on it, and produce a `.cif` file in the same location. Notice you do need pybabel, rdkit, etc to work. Details can be found in the `requirements.txt`.

### Detail of the Data Structure
This section below explains detail for each record following the same structure as saved in the .json files. 
#### struct_id
Material reference name, consistent with Cambridge Structural Database (CSD) reference name. 
#### geometry
Both the relaxed crystal structure and single molecule structure are saved under this key.  
__relaxed\_crystal__   
A dictionary format of pymatgen Structure.  
__molecule__  
A dictionary format of pymatgen Molecule.  
#### dft 
This section stores all the non-GWBSE values.  

```
bandgap: the crystal electronic band gap
Et: the crystal triplet formation energy, calculated by the total energy difference between the ground-state and triplet-state crystal
DF: the crystal DFT estimate for the SF driving force, calculated by taking the difference between bandgap and twice Et
VBdisp: the valence band dispersion, calculated by the energy range of the HOMO-derived band
CBdisp: the conduction band dispersion, calculated by the energy range of the LUMO-derived band
hab: the transfer integral, calculated with fragment orbital DFT
gap_s: the single molecule gap, calculated by the energy difference between highest occupied molecular orbital (HOMO) and lowest unoccupied molecular orbital (LUMO)
Et_s: the single molecule triplet formation energy, calculated by the total energy difference between the ground-state and triplet-state molecule
DF_s: the single molecule DFT estimate for the SF driving force, calculated by taking the difference between gap_s and twice Et_s
IP_s: the single molecule ionization potential, calculated by the total energy difference between a cation and neutral molecule
EA_s: the single molecule electron affinity, calculated by the total energy difference between an anion and neutral molecule
polarisation: the trace of the polarization tensor for a single molecule, calculated with DFT using the PBE functional and many-body dispersion (MBD) method (PBE+MBD)
apc: the number of atoms in the crystal unit cell
density: the crystal density in amu Å^−3
epsilon: the dielectric constant calculated with PBE+MBD
weight_s: the molecular weight in atomic mass unit (amu)
```
#### gwbse
__Es__  
The optical gap for singlet-state exciton   
__Et__  
The optical gap for triplet-state exciton  
__DF__  
The singlet fission driving force calculated by subtract twice Et from Es  
__absorption__  
Absorption has three directions `a`, `b`, and `c`, saved under the path `data["gwbse"]["absorption"]["a"]` `data["gwbse"]["absorption"]["b"]` `data["gwbse"]["absorption"]["c"]`, respectively. The four columns are, taken from BerkeleyGW output:   

```
 # Column 1: omega  
 # Column 2: eps2(omega)  
 # Column 3: eps1(omega)  
 # Column 4: JDOS(omega)  
```
The solar spectrum is also attached for comparison, under the path `data["gwbse"]["absorption"]["solar"]`  
__bandstructure__  
The direct absorption from BerkeleyGW is named as `bandstructure.dat` with first two lines commented as

```
# spin      band      kx      ky      kz      E(MF)      E(QP)      Delta E
#                  (Cartesian coordinates)     (eV)       (eV)         (eV)
   
```  
The colomn kx, ky, kz indicates the coordinates of each k points, the quasi-partical energy is given in the column E(QP). These four columns are used to plot the bandstructure. Bandstructure is saved under `data['gwbse']['bandstructure']['val']`. The high symmetry points used to build the k path are saved in `data['gwbse']['bandstructure']['kpoints']`. Both values are saved as list.    