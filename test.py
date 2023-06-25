# import json 
# import os
# from pymatgen.analysis.structure_matcher import StructureMatcher
# from ase.io import read, write
# from pymatgen.core import Structure
# from pymatgen.io.ase import AseAtomsAdaptor
# import matplotlib.pyplot as plt

# relaxed_struct_path = "/Users/alfred/Desktop/work/PAH101/PAH101Plot/data"
# unrelaxed_struct_path = "/Users/alfred/Desktop/work/PAH101/data/structs/crystal_unrelaxed"

# # load structures
# relaxed_struct = {}
# unrelaxed_struct = {}

# for name in os.listdir(relaxed_struct_path):
#     with open(os.path.join(relaxed_struct_path, name)) as f:
#         data = json.load(f)
#     relaxed_struct[name.split(".")[0]] = Structure.from_dict(data["geometry"]["relaxed_crystal"])

# for k, v in relaxed_struct.items():
#     relaxed_struct[k].remove_species("H")

# for name in os.listdir(unrelaxed_struct_path):
#     tmp_struct = read(os.path.join(unrelaxed_struct_path, name))
#     unrelaxed_struct[name.split(".")[0]] = AseAtomsAdaptor.get_structure(tmp_struct)

# for k, v in unrelaxed_struct.items():
#     unrelaxed_struct[k].remove_species("H")

# rmsd = []
# max_dist = []
# name_list = []

# # import code
# # code.interact(local=locals())

# for name in relaxed_struct.keys():
#     try:
#         print(name)
#         v1, v2 = StructureMatcher().get_rms_dist(struct1=relaxed_struct[name], struct2=unrelaxed_struct[name])
#         rmsd.append(v1)
#         max_dist.append(v2)
#         name_list.append(name)
#     except:
#         print("error")

# print("number of struct passed the test", len(rmsd))
# import code
# code.interact(local=locals())

# fig = plt.figure(num=None,figsize=(7, 5), dpi=80,facecolor='w', edgecolor='k', frameon=True)
# plt.hist(rmsd, 10, density = 1, color ='blue', alpha = 0.7)
# plt.axis([-0.01, 0.08, 0, 50])
# plt.xlabel('Root-Mean Squared Distance', fontsize=16)
# plt.ylabel('Count', fontsize=16)
# fig.savefig("rmsd.png", dpi=100, bbox_inches='tight')

# fig = plt.figure(num=None,figsize=(7, 5), dpi=80,facecolor='w', edgecolor='k', frameon=True)
# plt.hist(max_dist, 10, density = 1, color ='red', alpha = 0.7)
# plt.axis([-0.01, 0.35, 0, 20])
# plt.xlabel('Maximum Pair Distance', fontsize=16)
# plt.ylabel('Count', fontsize=16)
# fig.savefig("max_dist.png", dpi=100, bbox_inches='tight')

import os
import json
import numpy as np

relaxed_struct_path = "/Users/alfred/Desktop/work/PAH101/PAH101Plot/data"

data = {}
for name in os.listdir(relaxed_struct_path):
    with open(os.path.join(relaxed_struct_path, name)) as f:
        d = json.load(f)
    data[name.split(".")[0]] = d

target = "FLUANT02"
targetpath = "FLUANT02"

abs_a = np.loadtxt(os.path.join(targetpath, "a-absorption_eh.dat"))
abs_b = np.loadtxt(os.path.join(targetpath, "b-absorption_eh.dat"))
abs_c = np.loadtxt(os.path.join(targetpath, "c-absorption_eh.dat"))

data[target]["gwbse"]["absorption"]["a"] = abs_a.tolist()
data[target]["gwbse"]["absorption"]["b"] = abs_b.tolist()
data[target]["gwbse"]["absorption"]["c"] = abs_c.tolist()

banddata = np.loadtxt(os.path.join(targetpath, "bandstructure.dat"))
with open(os.path.join(targetpath, "kpoints")) as f:
    kpoints = f.readlines()
for i, n in enumerate(kpoints):
    # kpoints[i] = n[:-1]
    kpoints[i] = n.replace("\n", "")
    kpoints[i] = kpoints[i].replace("\t", "  ")

data[target]["gwbse"]["bandstructure"]["val"] = banddata.tolist()
data[target]["gwbse"]["bandstructure"]["kpoints"] = kpoints

print(data[target]["gwbse"]["bandstructure"]["kpoints"])

with open(target+"_new.json", "w") as f:
    json.dump(data[target], f, indent=4)










