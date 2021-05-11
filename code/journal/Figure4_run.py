from msmpb_util import *
import pandas as pd
import json
import matplotlib.pyplot as plt

experiment = {'data': {'name': 'Box44', 'N': 65850, 'A': 88},
       'trial': [ 
            {'nu':4, 'M': 90, 'a0': 9, 'msmout': {}, 'pmeout': {}}, 
            {'nu':4, 'M': 84, 'a0': 9, 'msmout': {}, 'pmeout': {}},
            {'nu':4, 'M': 80, 'a0': 9, 'msmout': {}, 'pmeout': {}},
            {'nu':4, 'M': 72, 'a0': 9, 'msmout': {}, 'pmeout': {}},
            {'nu':4, 'M': 70, 'a0': 9, 'msmout': {}, 'pmeout': {}},  
            {'nu':6, 'M': 60, 'a0': 9, 'msmout': {}, 'pmeout': {}}, 
            {'nu':6, 'M': 56, 'a0': 9, 'msmout': {}, 'pmeout': {}},
            {'nu':6, 'M': 54, 'a0': 9, 'msmout': {}, 'pmeout': {}},   
            {'nu':6, 'M': 50, 'a0': 9, 'msmout': {}, 'pmeout': {}}, 
            {'nu':6, 'M': 48, 'a0': 9, 'msmout': {}, 'pmeout': {}},                        
            {'nu':4, 'M': 54, 'a0':12, 'msmout': {}, 'pmeout': {}}, 
            {'nu':4, 'M': 50, 'a0':12, 'msmout': {}, 'pmeout': {}},   
            {'nu':4, 'M': 48, 'a0':12, 'msmout': {}, 'pmeout': {}}, 
            {'nu':4, 'M': 42, 'a0':12, 'msmout': {}, 'pmeout': {}}, 
            {'nu':4, 'M': 40, 'a0':12, 'msmout': {}, 'pmeout': {}},   
            {'nu':6, 'M': 42, 'a0':12, 'msmout': {}, 'pmeout': {}},  
            {'nu':6, 'M': 40, 'a0':12, 'msmout': {}, 'pmeout': {}}, 
            {'nu':6, 'M': 36, 'a0':12, 'msmout': {}, 'pmeout': {}},   
            {'nu':6, 'M': 32, 'a0':12, 'msmout': {}, 'pmeout': {}}, 
            {'nu':6, 'M': 30, 'a0':12, 'msmout': {}, 'pmeout': {}} ]}

for t in experiment['trial']:
    t['msmout'] = runMSM(experiment['data'],t['M'],t['a0'],t['nu'])
    t['pmeout'] = runPME(experiment['data'],t['M'],t['a0'],t['nu'])
    
with open('../../results/Figure4.json', 'w') as outfile:
    json.dump(experiment, outfile)
