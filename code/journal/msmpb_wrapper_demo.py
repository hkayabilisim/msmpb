from msmpb_wrapper import *  

data = {'name':'spceN300' ,'N':300  ,'A':20}

msm = runMSM(data, M=75, a0=4, nu=6, L=1)
pme = runPME(data, M=75, a0=4, nu=6)

print("Number of atoms %d " % msm['N'])
print("Energy (MSM) %f " % msm['utotal'])
print("Energy (PME) %f " % pme['utotal'])

