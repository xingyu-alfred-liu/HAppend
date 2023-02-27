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

### Data Detail
__Bandstructure__

__Absorption__  
Absorption has three directions `a`, `b`, and `c`, saved under the corresponding names. The four columns are, taken from BerkeleyGW output:   

```
 # Column 1: omega  
 # Column 2: eps2(omega)  
 # Column 3: eps1(omega)  
 # Column 4: JDOS(omega)  
```
The solar spectrum is also attached for comparison, under the path `data["gwbse"]["absorption"]["solar"]`.