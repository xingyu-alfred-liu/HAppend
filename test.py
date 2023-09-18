from pymatgen.core import Structure
from pymatgen.core import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
import json
import os
from ase.io import read, write

def determine_grid(x):
    if x <= 5:
        return 8
    if 5 < x <= 10:
        return 4
    if 10 < x <=20:
        return 2
    if x > 20:
        return 1

data_dict = {}
for name in os.listdir("data"):
    with open(os.path.join("data", name), "r") as f:
        data = json.load(f)
    m = Molecule.from_file(os.path.join("es", name.split(".")[0], "geo.xyz"))
    data["geometry"]["relaxed_molecule"] = m.to(fmt="json")
    with open(os.path.join("data_new", name), "w") as f:
        json.dump(data, f, indent=4)