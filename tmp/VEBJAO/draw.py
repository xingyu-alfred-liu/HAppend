import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import BSpline, make_interp_spline

# set the path & index from the command line
if __name__ == '__main__':
    a_inp = '4-a-absorption_eh.dat'  # path to the absorption_eh.dat file along a direction
    b_inp = '4-b-absorption_eh.dat'  # path to the absorption_eh.dat file along b direction
    c_inp = '4-c-absorption_eh.dat'  # path to the absorption_eh.dat file along c direction
    opt_inp = '/Users/siyugao/Desktop/research_files/singlet_fission/crystal_absorption/allsmrtetr.txt' # path to solar spectrum file
    # output = sys.argv[4] # path to output file
st_name = 'VEBJAO'
opt_gap = 2.05

# read omega as X axis, eps2(omega) as Y axis
def np_read(file): 
    X = []
    Y = []
    myfile = open(file)
    new = myfile.readlines()[4:]
    #  print(new)
    np.array(new)
    for line in new:
        X.append(float(line.split()[0])) 
        Y.append(float(line.split()[1]))
        # print(Y)
    Y_max = max(Y)
        # print(Y_max)
    Z = [i*500/Y_max for i in Y]
    myfile.close()
    return X, Z

# read solar spectrum
def sl_read(file): 
    X = []
    Y = []
    myfile = open(file)
    new = myfile.readlines()[1:2002]
    #  print(new)
    np.array(new)
    for line in new:
        X.append(1240.0/float(line.split()[0])) 
        Y.append(float(line.split()[1])* float(line.split()[0])/(1240.0/float(line.split()[0])))
    myfile.close()
    return np.flipud(X), np.flipud(Y)


(X, Y_a) = np_read(a_inp)
(X, Y_b) = np_read(b_inp)
(X, Y_c) = np_read(c_inp)
Y = np.add(np.add(Y_a, Y_b), Y_c)
(X_sl, Y_sl) = sl_read(opt_inp)

X_fill = []
Y_fill = []
for i, j in zip(X, Y):
    if i >opt_gap:
        X_fill.append(i)
        Y_fill.append(j)
#  [i in X if i > opt_gap ]

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(X, Y_a, '--k', label='GW/BSE-a')
ax1.plot(X, Y_b, '--m', label='GW/BSE-b')
ax1.plot(X, Y_c, '--y', label='GW/BSE-c')
ax1.plot(X, Y, '-b', label='Total')
ax1.fill_between(X_fill, Y_fill, 0)
# plt.plot(X_sl, Y_sl, '-c', label='Solar Spectrum')
ax1.axvline(x=opt_gap, linewidth=2, color='r',label='Optical Gap')
# x_new = np.linspace(1.5, 2.7, 300)
# a_BSpline = make_interp_spline(X_exp, Y_exp)
# y_new = a_BSpline(x_new)
# plt.plot(x_new, y_new, '-m', label='Experiment')
ax1.legend(loc='best',prop={'size': 14})
ax1.set_xlabel('Energy (eV)',fontsize=18)
ax1.yaxis.label.set_color('b')
ax1.set_ylabel('Intensity (arb.u.)',fontsize=18)
for tick in ax1.xaxis.get_majorticklabels():  # example for xaxis
    tick.set_fontsize(14) 
# ax1.set_xticks(fontsize=14)
ax1.set_yticks([])
plt.tick_params(axis='both', which='major', labelsize=14)
ax2.plot(X_sl, Y_sl, '-c', label='Solar Spectrum')
ax2.set_ylabel('Spectral Power (W $m^{-2}$ $eV^{-1}$)',fontsize=18)
ax2.yaxis.label.set_color('c')
ax2.tick_params(axis='y', colors='c')
# ax2.set_yticks(fontsize=14)
ax1.axis([1.25,8,0,800])
ax2.axis([1.25,8,0,800])
fig.set_size_inches(6.8, 5.1)
# plt.xlim(0,10)
plt.xlim(1.5,4.5)
plt.subplots_adjust(left=0.08, bottom=0.14, right=0.86, top=0.91)
plt.title('Absorption Spectrum %s' % st_name,fontsize=20)
# plt.show()
plt.savefig('Absorption Spectrum %s' % st_name)