import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json

experiments = [{"L": 1, "nu": 4, "M": [45, 50, 60, 64, 70, 72, 75, 80, 90]},
               {"L": 1, "nu": 6, "M": [45, 50, 60, 64, 70, 72, 75, 80, 90]},
               {"L": 2, "nu": 4, "M": [32, 36, 40, 42, 48]},
               {"L": 2, "nu": 6, "M": [32, 36, 40, 42, 48]},
               {"L": 3, "nu": 4, "M": [28, 32, 36, 40, 48]},
               {"L": 3, "nu": 6, "M": [28, 32, 36, 40, 48]}]
           
href={4:0.83,6:1}
aref={4:7   ,6:6}

print("%1s %2s %10s %4s %6s %8s %8s"%("L","nu","M(finest)","h","a0","time","error"))
for experiment in experiments:
    L = experiment["L"]
    nu = experiment["nu"]
    Mrange = experiment["M"]
    experiment["time"] = []
    experiment["error"] = []
    for M in Mrange:
        A = 88 # Box44 data
        h1 = A/M
        a1_0 = aref[nu]*((h1/href[nu])**(1-1/nu))
        a0 = a1_0 #* ((4/3)*(1-2**(-L)))**(1/nu)
        times = np.zeros((5,1))
        for run in range(5):
            out = subprocess.run(['./msmpb','../data/Box44','-L',str(L),'--a0',str(a0),'--nu',str(nu),'-M',str(M/(2**(L-1)))]\
                     ,stdout=subprocess.PIPE)
            msmout = json.loads(out.stdout)
            times[run]=msmout['time_total']
            error=msmout['deltaF/Fref']
        time = np.mean(times)
        experiment["time"].append(time)
        experiment["error"].append(error)
        print("%1d %2d %10d %4.2f %6.3f %8.3f %8.3e"%(L,  nu,  M, h1,   a0,   time,error))
