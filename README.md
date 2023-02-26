# PAH101Plot
Code to reproduce the plots in paper. In addition, the add Hydrogen program is attahced. 

### Data Detail
#### Bandstructure

#### Absorption
Absorption has three directions `a`, `b`, and `c`, saved under the corresponding names. The four columns are, taken from BerkeleyGW output:   

```
 # Column 1: omega  
 # Column 2: eps2(omega)  
 # Column 3: eps1(omega)  
 # Column 4: JDOS(omega)  
```
The solar spectrum is also attached for comparison, under the path `data["gwbse"]["absorption"]["solar"]`.