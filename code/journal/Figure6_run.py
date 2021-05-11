from msmpb_wrapper import *
import pandas as pd
import json
import matplotlib.pyplot as plt

experiment = {'data': {'name': 'Box44', 'N': 65850, 'A': 88},
       'trial': [ 
            {'L':1,'nu':4, 'M': 90, 'a0': 9, 'msmout': []},
            {'L':1,'nu':4, 'M': 84, 'a0': 9, 'msmout': []},
            {'L':1,'nu':4, 'M': 80, 'a0': 9, 'msmout': []},
            {'L':1,'nu':4, 'M': 72, 'a0': 9, 'msmout': []},
            {'L':1,'nu':4, 'M': 70, 'a0': 9, 'msmout': []}, 
            {'L':1,'nu':6, 'M': 60, 'a0': 9, 'msmout': []},
            {'L':1,'nu':6, 'M': 56, 'a0': 9, 'msmout': []},
            {'L':1,'nu':6, 'M': 54, 'a0': 9, 'msmout': []},  
            {'L':1,'nu':6, 'M': 50, 'a0': 9, 'msmout': []},
            {'L':1,'nu':6, 'M': 48, 'a0': 9, 'msmout': []},
            {'L':1,'nu':4, 'M': 54, 'a0':12, 'msmout': []},
            {'L':1,'nu':4, 'M': 50, 'a0':12, 'msmout': []},  
            {'L':1,'nu':4, 'M': 48, 'a0':12, 'msmout': []},
            {'L':1,'nu':4, 'M': 42, 'a0':12, 'msmout': []},
            {'L':1,'nu':4, 'M': 40, 'a0':12, 'msmout': []},  
            {'L':1,'nu':6, 'M': 42, 'a0':12, 'msmout': []}, 
            {'L':1,'nu':6, 'M': 40, 'a0':12, 'msmout': []},
            {'L':1,'nu':6, 'M': 36, 'a0':12, 'msmout': []},  
            {'L':1,'nu':6, 'M': 32, 'a0':12, 'msmout': []},
            {'L':1,'nu':6, 'M': 30, 'a0':12, 'msmout': []},
            {'L':2,'nu':6, 'M': 60, 'a0': 9, 'msmout': []},
            {'L':2,'nu':6, 'M': 56, 'a0': 9, 'msmout': []},
            {'L':2,'nu':6, 'M': 54, 'a0': 9, 'msmout': []},  
            {'L':2,'nu':6, 'M': 50, 'a0': 9, 'msmout': []},
            {'L':2,'nu':6, 'M': 48, 'a0': 9, 'msmout': []},
            {'L':2,'nu':4, 'M': 54, 'a0':12, 'msmout': []},
            {'L':2,'nu':4, 'M': 50, 'a0':12, 'msmout': []},  
            {'L':2,'nu':4, 'M': 48, 'a0':12, 'msmout': []},
            {'L':2,'nu':4, 'M': 42, 'a0':12, 'msmout': []},
            {'L':2,'nu':4, 'M': 40, 'a0':12, 'msmout': []},  
            {'L':2,'nu':6, 'M': 42, 'a0':12, 'msmout': []}, 
            {'L':2,'nu':6, 'M': 40, 'a0':12, 'msmout': []},
            {'L':2,'nu':6, 'M': 36, 'a0':12, 'msmout': []},  
            {'L':2,'nu':6, 'M': 32, 'a0':12, 'msmout': []},
            {'L':2,'nu':6, 'M': 30, 'a0':12, 'msmout': []},
            {'L':3,'nu':6, 'M': 60, 'a0': 9, 'msmout': []},
            {'L':3,'nu':6, 'M': 56, 'a0': 9, 'msmout': []},
            {'L':3,'nu':6, 'M': 48, 'a0': 9, 'msmout': []},
            {'L':3,'nu':4, 'M': 48, 'a0':12, 'msmout': []},
            {'L':3,'nu':4, 'M': 40, 'a0':12, 'msmout': []},  
            {'L':3,'nu':6, 'M': 36, 'a0':12, 'msmout': []},  
            {'L':3,'nu':6, 'M': 32, 'a0':12, 'msmout': []}]}

for t in experiment['trial']:
    topLevelM = t['M']/(2**(t['L']-1))
    t['msmout'] = runMSM(experiment['data'],topLevelM,t['a0'],t['nu'],L=t['L'])
    print(t['msmout'])
    if t['L'] != t['msmout']['NumberOfLevels'] or t['nu'] != t['msmout']['nu'] or t['a0'] != t['msmout']['cutoff'] or topLevelM != t['msmout']['TopLevelMx']:
       print("\nError\n")
    
with open('../../results/Figure6.json', 'w') as outfile:
    json.dump(experiment, outfile)

