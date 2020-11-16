import numpy as np
import matplotlib.pyplot as plt
from msmpb_util import *
import subprocess
import json


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
displayResult(training,fudgefactors)
print("Testing -- (fudgefactors applied!)")
displayResult(testing,fudgefactors)


plt.figure(figsize=(40,20))
plt.rcParams.update({'font.size': 22})
ax1 = plt.subplot(1,2,1)
plotEstimatedVersusTrue(ax1,training,fudgefactors)
ax2 = plt.subplot(1,2,2)
plotEstimatedVersusTrue(ax2,testing,fudgefactors)
plt.savefig('../results/Figure2-nu.pdf')  

plt.figure(figsize=(40,20))
plt.rcParams.update({'font.size': 22})
ax1 = plt.subplot(1,2,1)
plotEstimatedVersusTrueWRTa0(ax1,training,fudgefactors)
ax2 = plt.subplot(1,2,2)
plotEstimatedVersusTrueWRTa0(ax2,testing,fudgefactors)
plt.savefig('../results/Figure2-a0.pdf')  

plt.figure(figsize=(40,20))
plt.rcParams.update({'font.size': 22})
ax1 = plt.subplot(1,2,1)
plotEstimatedVersusTrueWRTM(ax1,training,fudgefactors)
ax2 = plt.subplot(1,2,2)
plotEstimatedVersusTrueWRTM(ax2,testing,fudgefactors)
plt.savefig('../results/Figure2-M.pdf')  

plt.figure(figsize=(40,20))
plt.rcParams.update({'font.size': 22})
ax1 = plt.subplot(1,2,1)
plotEstimatedVersusTrueWRTN(ax1,training,fudgefactors)
ax2 = plt.subplot(1,2,2)
plotEstimatedVersusTrueWRTN(ax2,testing,fudgefactors)
plt.savefig('../results/Figure2-N.pdf')  
