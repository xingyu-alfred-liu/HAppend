import json
import numpy as np
import os

for fname in os.listdir("bandstructure_dat"):
    if "kpoints" in os.listdir("bandstructure_dat/%s" % (fname)):
        with open('data/%s.json' % (fname), "r") as f:
            data = json.load(f)
        # abs_a = np.loadtxt("bandstructure_dat/%s/4-a-absorption_eh.dat" % (fname))
        # abs_b = np.loadtxt("bandstructure_dat/%s/4-b-absorption_eh.dat" % (fname))
        # abs_c = np.loadtxt("bandstructure_dat/%s/4-c-absorption_eh.dat" % (fname))
        bands_data = np.loadtxt("bandstructure_dat/%s/bandstructure.dat" % (fname))
        with open("bandstructure_dat/%s/kpoints" % (fname),"r") as f:
            kpoints = f.readlines()
        for i, v in enumerate(kpoints):
            kpoints[i] = v.replace("\n", "")
        for i, v in enumerate(kpoints):
            kpoints[i] = v.replace("\t", "  ")
        # data["gwbse"]["absorption"]["a"] = abs_a.tolist()
        # data["gwbse"]["absorption"]["b"] = abs_b.tolist()
        # data["gwbse"]["absorption"]["c"] = abs_c.tolist()
        data["gwbse"]["bandstructure"]["val"] = bands_data.tolist()
        data["gwbse"]["bandstructure"]["kpoints"] = kpoints
        with open("bandstructure_dat/%s.json" % (fname), "w") as f:
            json.dump(data, f, indent=4)