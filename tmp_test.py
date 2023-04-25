import json
import numpy as np
fname = "VEBKAP"
with open('data/%s.json' % (fname), "r") as f:
    data = json.load(f)
abs_a = np.loadtxt("tmp/%s/4-a-absorption_eh.dat" % (fname))
abs_b = np.loadtxt("tmp/%s/4-b-absorption_eh.dat" % (fname))
abs_c = np.loadtxt("tmp/%s/4-c-absorption_eh.dat" % (fname))
# bands_data = np.loadtxt("tmp/%s/bandstructure.dat" % (fname))
# with open("tmp/%s/kpoints" % (fname),"r") as f:
#     kpoints = f.readlines()
# for i, v in enumerate(kpoints):
#     kpoints[i] = v.replace("\n", "")
# for i, v in enumerate(kpoints):
#     kpoints[i] = v.replace("\t", "  ")
data["gwbse"]["absorption"]["a"] = abs_a.tolist()
data["gwbse"]["absorption"]["b"] = abs_b.tolist()
data["gwbse"]["absorption"]["c"] = abs_c.tolist()
# data["gwbse"]["bandstructure"]["val"] = bands_data.tolist()
# data["gwbse"]["bandstructure"]["kpoints"] = kpoints
with open("data/%s_new.json" % (fname), "w") as f:
    json.dump(data, f, indent=4)