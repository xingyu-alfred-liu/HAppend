import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
import json

with open("BOXGAW.json", "r") as f:
  boxgaw = json.load(f)
boxgaw_geo = read("../BOXGAW/BOXGAW.cif")
boxgaw_geo = AseAtomsAdaptor.get_structure(boxgaw_geo)
boxgaw["geometry"]["relaxed_crystal"] = boxgaw_geo.as_dict()
# bandstructure
boxgaw_band = np.loadtxt("../BOXGAW/bandstructure.dat")
with open ("../BOXGAW/kpoints", "r") as f:
  boxgaw_kpoints = f.readlines()
for i, d in enumerate(boxgaw_kpoints):
  boxgaw_kpoints[i] = d[:-1]
boxgaw["gwbse"]["bandstructure"]["kpoints"] = boxgaw_kpoints
boxgaw["gwbse"]["bandstructure"]["val"] = boxgaw_band.tolist()
boxgaw_a = np.loadtxt("../BOXGAW/a-absorption_eh.dat")
boxgaw_b = np.loadtxt("../BOXGAW/b-absorption_eh.dat")
boxgaw_c = np.loadtxt("../BOXGAW/c-absorption_eh.dat")
boxgaw["gwbse"]["absorption"]["a"] = boxgaw_a.tolist()
boxgaw["gwbse"]["absorption"]["b"] = boxgaw_b.tolist()
boxgaw["gwbse"]["absorption"]["c"] = boxgaw_c.tolist()
with open("BOXGAW_new.json", "w") as f:
  json.dump(boxgaw, f, indent=4)

with open("CENYEZ.json", "r") as f:
  cenyez = json.load(f)
cenyez_geo = read("../CENYEZ/CENYEZ.cif")
cenyez_geo = AseAtomsAdaptor.get_structure(cenyez_geo)
cenyez["geometry"]["relaxed_crystal"] = cenyez_geo.as_dict()
# bandstructure
cenyez_band = np.loadtxt("../CENYEZ/bandstructure.dat")
with open ("../CENYEZ/kpoints", "r") as f:
  cenyez_kpoints = f.readlines()
for i, d in enumerate(cenyez_kpoints):
  cenyez_kpoints[i] = d[:-1]
cenyez["gwbse"]["bandstructure"]["kpoints"] = cenyez_kpoints
cenyez["gwbse"]["bandstructure"]["val"] = cenyez_band.tolist()
# absorption
cenyez_a = np.loadtxt("../CENYEZ/a-absorption_eh.dat")
cenyez_b = np.loadtxt("../CENYEZ/b-absorption_eh.dat")
cenyez_c = np.loadtxt("../CENYEZ/c-absorption_eh.dat")
cenyez["gwbse"]["absorption"]["a"] = cenyez_a.tolist()
cenyez["gwbse"]["absorption"]["b"] = cenyez_b.tolist()
cenyez["gwbse"]["absorption"]["c"] = cenyez_c.tolist()
with open("CENYEZ_new.json", "w") as f:
  json.dump(cenyez, f, indent=4)