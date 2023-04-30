import numpy as np
import code

data = np.loadtxt('bandstructure.dat')

band_val = data[0][1]
for idx in range(data.shape[0]):
    if data[idx][1] != band_val:
        break
data_test = data[:idx]

kpoint_list = [data_test[0][2:5]]
kpoint_idx = [0]
dist = np.linalg.norm(data_test[1][2:5]-data_test[0][2:5])
kpoint_name = ["Gamma"]
idx = 0

while idx+1 < data_test.shape[0]:
    dist_compare = np.linalg.norm(data_test[idx+1][2:5]-data_test[idx][2:5])
    if dist_compare == 0:
        kpoint_idx.append(idx)
        kpoint_list.append(data_test[idx][2:5])
        kpoint_name.append("Gamma")
        
    if  abs(dist_compare - dist) > 1e-4:
        kpoint_idx.append(idx+1)
        kpoint_list.append(data_test[idx+1][2:5])
        kpoint_name.append("Gamma")
        if idx+2 < data_test.shape[0]:
            dist = np.linalg.norm(data_test[idx+2][2:5]-data_test[idx+1][2:5])
            idx = idx+2
    else:
        idx += 1 

for i in range(len(kpoint_idx)-1):
    print("%s %s %s " % (str(kpoint_list[i][0]), str(kpoint_list[i][1]), str(kpoint_list[i][2])), "##%s" % kpoint_name[i])
    print(kpoint_idx[i+1]-kpoint_idx[i-1]-1)
print("%s %s %s " % (str(kpoint_list[-1][0]), str(kpoint_list[-1][1]), str(kpoint_list[-1][2])), "##%s" % kpoint_name[-1])


with open("kpoints", "w") as f:
    for i in range(len(kpoint_idx)-1):
        f.write(str(kpoint_list[i][0])+"  "+str(kpoint_list[i][1])+"  "+str(kpoint_list[i][2])+"  "+"##%s" % kpoint_name[i]+"\n")
        f.write(str(kpoint_idx[i+1]-kpoint_idx[i]-1)+"\n")
    f.write(str(kpoint_list[-1][0])+"  "+str(kpoint_list[-1][1])+"  "+str(kpoint_list[-1][2])+"  "+"##%s" % kpoint_name[i]+"\n")

code.interact(local=locals())