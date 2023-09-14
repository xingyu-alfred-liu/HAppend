import os
import numpy as np

data = {}
for name in os.listdir("es"):
    print(name)
    with open(os.path.join("es", name, "aims.out")) as f:
        p = f.readlines()
    data_record = []
    start_record = False
    for l in p:
        if "State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]" in l:
            start_record = True
        if "Highest occupied state (VBM) at" in l:
            start_record = False
        if start_record is True:
            data_record.append(l[:-1].split())
    new_data = []
    for i, l in enumerate(data_record):
        if len(l) == 4:
            new_data.append(l)
    new_data_record = np.array(new_data)
    
    final_record = []
    for i in range(len(new_data_record)-1, -1, -1):
        final_record.append(new_data_record[i])
        if new_data_record[i][0] == '1':
            break
    final_record = np.array(final_record[::-1])
    print(final_record.shape)
    data[name] = final_record.astype(float)

import json
for k, v in data.items():
    print(k)
    with open("../data/%s.json" % (k), "r") as f:
        tmp_data = json.load(f)
    tmp_data["dft"]["eigenvalues"] = v.tolist()

    with open("%s_new.json" % (k), "w") as f:
        json.dump(tmp_data, f, indent=4)