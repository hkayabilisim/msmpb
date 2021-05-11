from msmpb_wrapper import *
import pandas as pd
import json
import matplotlib.pyplot as plt

experiment = {'data': {'name': 'Box44', 'N': 65850, 'A': 88},
       'trial': [ 
            {'L':2,'nu':6, 'M': 36, 'a0':12, 'msmout': []}]}

for t in experiment['trial']:
    topLevelM = t['M']/(2**(t['L']-1))
    t['msmout'] = runMSM(experiment['data'],topLevelM,t['a0'],t['nu'],L=t['L'])
    print(t['msmout'])
    if t['L'] != t['msmout']['NumberOfLevels'] or t['nu'] != t['msmout']['nu'] or t['a0'] != t['msmout']['cutoff'] or topLevelM != t['msmout']['TopLevelMx']:
       print("\nError-------\n")
    
with open('../../results/TableII.json', 'w') as outfile:
    json.dump(experiment, outfile)

