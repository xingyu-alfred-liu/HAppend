import numpy as np
import matplotlib.pyplot as plt
import json
with open("../data/TPHBEN01.json", "r") as f:
    data = json.load(f)
band = np.array(data["gwbse"]["bandstructure"]["val"])
xyz = band[:, 1:4]
val = band[:, 6].reshape(-1, 46)
fig = plt.figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', 
                edgecolor='k', frameon=True)
ax = fig.add_subplot(111)
for i in range(60):
    plt.plot(np.arange(46), val[i], 'r', linewidth=1.5)
    # plt.scatter(np.arange(60), val, s=6, marker='o', c='red')
    plt.savefig('test.png', bbox_inches='tight')