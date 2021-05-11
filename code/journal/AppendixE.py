from msmpb_wrapper import *

import sys

sys.stdout = open("../../results/AppendixE.txt", "w")

experiments = [{"L": 1, "nu": 4, "M": [80, 84, 90, 96,100]},
               {"L": 1, "nu": 6, "M": [80, 84, 90, 96,100]},
               {"L": 2, "nu": 4, "M": [40, 42, 48, 50, 54]},
               {"L": 2, "nu": 6, "M": [40, 42, 48, 50, 54]},
               {"L": 3, "nu": 4, "M": [36, 40, 48, 56, 60]},
               {"L": 3, "nu": 6, "M": [36, 40, 48, 56, 60]}]

href={4:0.83,6:1}
aref={4:7   ,6:6}

box44 = {'name':'Box44'    ,'N':65850,'A':88}

print("%1s %2s %10s %4s %6s %8s %8s"%("L","nu","M(finest)","h","a0","time","error"))
for experiment in experiments:
    L = experiment["L"]
    nu = experiment["nu"]
    Mrange = experiment["M"]
    experiment["time"] = []
    experiment["error"] = []
    for M in Mrange:
        A = box44['A']
        h1 = A/M
        a1_0 = aref[nu]*((h1/href[nu])**(1-1/nu))
        a0 = a1_0 
        times = np.zeros((5,1))
        for run in range(5):
            Mtoplevel = M//(2**(L-1))
            msmout = runMSM(box44,Mtoplevel,a0,nu,L)
            times[run]=msmout['time_total']
            error=msmout['deltaF/Fref']
        time = np.mean(times)
        experiment["time"].append(time)
        experiment["error"].append(error)
        print("%1d %2d %10d %4.2f %6.3f %8.3f %8.3e"%(L,  nu,  M, h1,   a0,   time,error))
        
        
for order in [4,6]:
    print("Order is %d"%order)
    for level in [1,2,3]:
        for experiment in experiments:
            if experiment['L'] == level and experiment['nu'] == order:               
                print("%20s MSM L=%d "%("Force Errorx1e4",level),end='')
                for error in experiment['error']:
                    print("%5.2f "%(error*1e4),end='')
                print("")    
    for level in [1,2,3]:
        for experiment in experiments:
            if experiment['L'] == level and experiment['nu'] == order:                      
                print("%20s MSM L=%d "%("Time",level),end='')
                for error in experiment['time']:
                    print("%5.2f "%error,end='')
                print("")          
              
                
    
sys.stdout.close()
