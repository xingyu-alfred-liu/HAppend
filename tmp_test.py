import json
import numpy as np
import os
import matplotlib.pyplot as plt

for fname in os.listdir("bandstructure_dat"):
    if "kpoints" in os.listdir("bandstructure_dat/%s" % (fname)):
        with open('data/%s.json' % (fname), "r") as f:
            data = json.load(f)
        abs_a = np.loadtxt("bandstructure_dat/%s/a-absorption_eh.dat" % (fname))
        abs_b = np.loadtxt("bandstructure_dat/%s/b-absorption_eh.dat" % (fname))
        abs_c = np.loadtxt("bandstructure_dat/%s/c-absorption_eh.dat" % (fname))
        bands_data = np.loadtxt("bandstructure_dat/%s/bandstructure.dat" % (fname))
        with open("bandstructure_dat/%s/kpoints" % (fname),"r") as f:
            kpoints = f.readlines()
        for i, v in enumerate(kpoints):
            kpoints[i] = v.replace("\n", "")
        for i, v in enumerate(kpoints):
            kpoints[i] = v.replace("\t", "  ")
        data["gwbse"]["absorption"]["a"] = abs_a.tolist()
        data["gwbse"]["absorption"]["b"] = abs_b.tolist()
        data["gwbse"]["absorption"]["c"] = abs_c.tolist()
        data["gwbse"]["bandstructure"]["val"] = bands_data.tolist()
        data["gwbse"]["bandstructure"]["kpoints"] = kpoints
        with open("bandstructure_dat/%s.json" % (fname), "w") as f:
            json.dump(data, f, indent=4)



# fig = plt.figure(num=None, figsize=(7, 6), dpi=80, facecolor='w', \
#                  edgecolor='k', frameon=True)

# x_list = [350,450,550,650]
# y_list = [2.380953,2.324623,2.288795,2.279896]

# plt.plot(x_list, y_list, '-', markersize=6, color="blue")
# plt.plot(x_list, y_list, '.', markersize=9, color="blue")

# plt.xlabel('Number of Empty States', fontsize=16)
# plt.ylabel('Quasi-Particle Energy Gap (eV)', fontsize=16)
# plt.axis([300, 700, 2.25, 2.4])

# plt.tick_params(
#     axis='x',            # changes apply to the x-axis
#     which='both',        # both major and minor ticks are affected
#     labelsize=16,
#     length = 4,
#     labelbottom=True)   # labels along the bottom edge are off
# # y ticks
# plt.tick_params(
#     axis='y',            # changes apply to the x-axis
#     which='both',        # both major and minor ticks are affected
#     length = 4, 
#     direction=None,
#     labelsize=16)

# fig.savefig('convergence_test.png', dpi=300, bbox_inches='tight')