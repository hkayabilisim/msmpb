import numpy as np
import matplotlib.pyplot as plt
from msmpb_wrapper import *
import subprocess
import json
import pandas as pd

data = [{'name':'spceN300' ,'N':300  ,'A':20},
        {'name':'spceN600' ,'N':600  ,'A':20},
        {'name':'spceN900' ,'N':900  ,'A':20},
        {'name':'spceN2250','N':2250 ,'A':30},
        {'name':'Box44'    ,'N':65850,'A':88}]

eta = 5.7

training = []
testing = []

print("Before calculating fudge factor")
print("%20s %5s %2s %5s %2s %8s %2s %10s %10s %s"%("data","N","A","nbar","M","a0","nu","true","est","|log(true)-log(est)|"))

for benchmark in data:
    N = benchmark['N']
    A = benchmark['A']
    M = makeFFTFriendly(calculateM(eta,N),decrease=False)
    Mnext = nextFFTFriendly(M)
    Mprev = prevFFTFriendly(M)
    for nbar in [12,14,16,18,20,22]:
        a0 = nbar * (A/M) / 2
        for nu in [4,6,8,10]:
            out = runMSM(benchmark,M,a0,nu)
            est = out["deltaFest/Fref"];
            true = out["deltaF/Fref"]
            diff = np.abs(np.log(est)-np.log(true))
            print("%20s %5d %2d %5.2f %2d %8.3f %2d %10.3e %10.3e %5.3f" \
              %(benchmark['name'],N,A,2 * a0 * M / A,M,a0,nu,true,est,diff))
            training.append(out)
            testing.append(runMSM(benchmark,Mprev,a0,nu))
            testing.append(runMSM(benchmark,Mnext,a0,nu))


fudgefactors = fudgeFactor(training)
print("Fudge factors (training)")
print(fudgefactors)

print("Training -- (fudgefactors applied!)")
displayResult(training,fudgefactors,"../../results/Figure2-training.txt")
print("Testing -- (fudgefactors applied!)")
displayResult(testing,fudgefactors,"../../results/Figure2-testing.txt")


trn = pd.read_csv('../../results/Figure2-training.txt',delim_whitespace=True)
tst = pd.read_csv('../../results/Figure2-testing.txt',delim_whitespace=True)

trntst = pd.concat([trn,tst])
data = [trn,trntst]


fig = plt.figure(figsize=(12, 4))
spec = fig.add_gridspec(ncols=2, nrows=1)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7)
for i in [0,1]:
    ax = fig.add_subplot(spec[0,i])
    subdata = data[i]
    ax.grid(b=True, which='major',linestyle='dotted',color='lightgrey')
    ax.plot([1e-8, 1e-2],[1e-8, 1e-2],'k',label=None)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("estimated")
    ax.set_ylabel("true")
    ax.set_xlim(1e-8, 1e-2)
    ax.set_ylim(1e-8, 1e-2)
    for nu in [4,6,8,10]:
        x = subdata[subdata['nu'] == nu]
        ax.scatter(x['est'],x['true'],label=r'$\nu='+str(nu)+'$',
                  edgecolor='k')
    ax.legend()
    #ax.axis('equal')

fig.savefig('../../results/Figure2.pdf', dpi=100)
 
