"""example to gather bandstructure and absorption
command: python plot_bandstructure_absorption.py -p ../data -id ABECAL -b True -a True -bp band.png -ap absorp.png
"""

import os
import json
import math

import numpy as np
import matplotlib.pyplot as plt

class PlotAgent:
    def __init__(self, root_path):
        self.data_dict = self.load_all(root_path)

    def load_all(self, root_path):
        """load all data into a dictionary
        Args:
            root_path (str): path to where all json are located
        Returns: 
            data_dict (Dict): dictionary of all data with name as key and 
            all content as value
        """
        data_dict = {}

        for name in os.listdir(root_path):
            if not name.endswith(".json"):
                continue

            with open(os.path.join(root_path, name), "r") as f:
                data = json.load(f)
                data_dict[name.split(".")[0]] = data

        return data_dict

    def plot_absorption(self, struct_id="ABECAL", savefig_path=None, ylim=20):
        data = self.data_dict[struct_id]["gwbse"]["absorption"]
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        absorb_a = np.array(data["a"])
        absorb_b = np.array(data["b"])
        absorb_c = np.array(data["c"])
        absorb_sum = absorb_a+absorb_b+absorb_c
        ax1.plot(absorb_a[:, 0], absorb_a[:, 1], '-b', label='Cal-a')
        ax1.plot(absorb_b[:, 0], absorb_b[:, 1], '-g', label='Cal-b')
        ax1.plot(absorb_c[:, 0], absorb_c[:, 1], '-m', label='Cal-c')
        ax1.plot(absorb_sum[:, 0], absorb_sum[:, 1], '-r', label='Cal-sum')
        # ax1.plot(X_exp_a, y_exp_a, '--b', label = 'Exp-a', linewidth=2)
        # ax1.plot(X_exp_b, y_exp_b, '--g', label = 'Exp-b', linewidth=2)
        
        ax1.legend(loc='upper right', prop={'size': 14})
        ax1.set_xlabel('Energy (eV)',fontsize=18)
        ax1.set_ylabel('Intensity (arb.u.)',fontsize=18)
        for tick in ax1.xaxis.get_majorticklabels():  # Adjust x ticks fontsize
            tick.set_fontsize(14) 
        ax1.set_yticks([])
        plt.tick_params(axis='both', which='major', labelsize=14)

        solar_val = np.array(data["solar"])
        ax2.plot(solar_val[:, 0], solar_val[:, 1], '-c', label='Solar Spectrum')
        ax2.set_ylabel('Spectral Power (W $m^{-2}$ $eV^{-1}$)',fontsize=18)
        ax2.yaxis.label.set_color('c')
        ax2.tick_params(axis='y', colors='c')
        # ax2.set_yticks(fontsize=14)
        ax1.axis([0.0, 10, 0, ylim])
        ax2.axis([0.0, 10, 0 ,800])
        fig.set_size_inches(6.8, 5.1)
        plt.xlim(0.0,10.0)
        plt.subplots_adjust(left=0.08, bottom=0.14, right=0.86, top=0.91)
        plt.title('Absorption Spectrum - %s' % (struct_id),fontsize=20)
        
        if savefig_path is not None:
            plt.savefig(savefig_path)
        
        return fig

    def plot_bandstructure(self, struct_id="ABECAL", savefig_path=None):
        data = self.data_dict[struct_id]["gwbse"]["bandstructure"]
        banddata = np.array(data["val"])
        kpoints = data["kpoints"]
        interList = []
        highsimPoint = []
        for i, val in enumerate(kpoints):
            line = val.split()
            if len(line) == 0:
                break
            elif len(line) == 1:
                interList.append(int(line[0]))
            else:
                name = line[-1]
                highsimPoint.append(name)
        print("The high symmetry points are:", highsimPoint)
        print("The invervals between each pair of points are:", interList)

        # now get the shift point index and the kpoint index
        shiftIndex = []
        kpointIndex = [0]
        counter = 0
        for i, val in enumerate(interList):
            if val != 0:
                counter += (val + 1)
                kpointIndex.append(counter)
            else:
                counter += 1
                shiftIndex.append(counter)
        print("The shift index are:", shiftIndex)
        print('The kpoint index are:', kpointIndex)
        kpointName = []
        i = 0
        while i < len(interList):
            if interList[i] != 0:
                kpointName.append(highsimPoint[i])
                i += 1
            else:
                kpointName.append(highsimPoint[i]+highsimPoint[i+1])
                i += 2
        kpointName.append(highsimPoint[-1])
        print("The kpoint names to be printed are:", kpointName)

        # get the number of points in one band
        numPoint = 1
        for i in range(len(banddata)):
            if banddata[i+1][1] != banddata[i][1]:
                break
            else:
                numPoint += 1
        # create a data array for band data
        QPband = np.zeros(shape=(numPoint, 1))
        # this is the list for special points

        for i in range(numPoint):
            if i == 0:
                QPband[i][0] == 0
            elif i in shiftIndex:
                QPband[i][0] = QPband[i-1][0]
            else:
                QPband[i][0]=math.sqrt((abs(banddata[i][2])-abs(banddata[i-1][2]))**2+
                    (abs(banddata[i][3])-abs(banddata[i-1][3]))**2 +
                    (abs(banddata[i][4])-abs(banddata[i-1][4]))**2) + QPband[i-1][0]
        # insert the band values, this is the QP band
        for i in range(int(len(banddata)/numPoint)):
            newcol = banddata[i*numPoint:(i+1)*numPoint, 6]
            QPband = np.append(QPband, newcol.reshape(numPoint, 1), 1)

        # Solve for band inversion by reordering the bands
        QP_band=QPband.T[1:]
        QP_band=np.sort(QP_band,axis=0)
        QPband=np.hstack((np.reshape(QPband[:,0],(-1,1)),QP_band.T))

        # the plot should be shift down shiftVal together
        row, col = QPband.shape
        shiftVal = max(QPband[:, int((col-1)/2)])

        # mark the band gap, indicate the max/min for homo/lumo
        lowGapIndex = np.argmax(QPband[:, int((col-1)/2)])
        highGapIndex = np.argmin(QPband[:, int((col-1)/2)+1])
        bandgapVal = max(QPband[:, int((col-1)/2)]) - min(QPband[:, int((col-1)/2)+1])

        print('bandgap is:', abs(bandgapVal), 'eV')
        # start the plot here
        fig = plt.figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', 
                        edgecolor='k', frameon=True)
        ax = fig.add_subplot(111)
        for i in range(len(QPband[0])-1):
            plt.plot(QPband[:,0], QPband[:,i+1]-shiftVal, 'r', linewidth=1.5)
            plt.scatter(QPband[:,0], QPband[:,i+1]-shiftVal, s=6, marker='o', c='red')

        # add the band gap plot
        plt.plot([QPband[lowGapIndex][0], QPband[highGapIndex][0]], \
                [QPband[lowGapIndex][int((col-1)/2)]-shiftVal, \
                QPband[highGapIndex][int((col-1)/2)+1]-shiftVal], '--b', linewidth=2)

        # add the band gap value
        plt.text(0.5*(QPband[lowGapIndex][0]+QPband[highGapIndex][0])+0.06, \
            0.5*(QPband[lowGapIndex][int((col-1)/2)]+QPband[highGapIndex][int((col-1)/2)+1])- 
            shiftVal+0.02, str(abs(bandgapVal))[:4]+' eV', weight='normal', 
            size='xx-large', color='blue')
        # add the high symmtry point
        plt.axis([0.0, QPband[numPoint-1][0], -3, 7])

        # manipulate the kpoint into $*$ format
        for i, name in enumerate(kpointName):
            print("name", name)
            if "Gamma" in name:
                kpointName[i] = "$" + "\\" + name[2:] + "$"
            else:
                kpointName[i] = "$" + name[2:] + "$"
        print(kpointName)
        # add high symmetry points
        for i in kpointIndex:
            plt.plot([QPband[i][0], QPband[i][0]], [-100, 100], 'k', linewidth=2)
        # add the fermi energy
        plt.plot([0.0, QPband[numPoint-1][0]], [0.0, 0.0], 'k--', linewidth=1.0)
        
        # add x y label
        plt.ylabel('$E-E_{F}$ (eV)', fontname="Arial", fontsize=20)
        nameVal = []
        for i in kpointIndex:
            nameVal += [QPband[i][0]]
        nameVal.sort()
        # x ticks
        plt.tick_params(
            axis='x',            # changes apply to the x-axis
            which='both',        # both major and minor ticks are affected
            bottom=False,        # ticks along the bottom edge are off
            top=False,           # ticks along the top edge are off
            labelbottom=True)   # labels along the bottom edge are off
        plt.xticks(nameVal, kpointName, fontname="Arial", fontsize=14)
        # y ticks
        plt.tick_params(
            axis='y',            # changes apply to the x-axis
            which='both',        # both major and minor ticks are affected
            bottom=False,        # ticks along the bottom edge are off
            top=False,           # ticks along the top edge are off
            labelbottom=True)    # labels along the bottom edge are off
        plt.yticks(fontname="Helvatica", fontsize=18)
        plt.title(struct_id, fontname='Helvatica', fontsize=22)

        if savefig_path is not None:
            fig.savefig(savefig_path, dpi=300, bbox_inches='tight')

        return fig

# if __name__ == "__main__":
#     main()