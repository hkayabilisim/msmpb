from msmpb_util import *
import pandas as pd
import json
import matplotlib.pyplot as plt

with open('../../results/TableII.json', 'r') as inpfile:
    experiment = json.load(inpfile)

items = {'time_stencil_atlevel1': 'Stencil at l=1',
         'time_stencil_atlevel2': 'Stencil at l=2',
         'time_nlist'           : 'Neighbor list',
         'time_partcl2partcl'   : 'Particle-to-particle',
         'time_anterpolation'   : 'Anterpolation',
         'time_interpolation'   : 'Interpolation',
         'time_grid2grid_atlevel1'  : 'Grid-to-grid at l=1',
         'time_restrictionfrom1to2'   : 'Restriction l=1->l=2',
         'time_prolongationfrom2to1'  : 'Prolongation l=2->l=1',
         'time_grid2grid_atlevel2'  : 'Grid-to-grid at l=2 (FFT)'}
print("%-30s " % "Phase",end='')
for t in experiment['trial']:
  L=t['msmout']['NumberOfLevels']
  TopM=t['msmout']['TopLevelMx']
  M = TopM * (2**(L-1))
  print("& %7s=%2s " % ('M',M),end='') 
print()
for item in items:
  print("%-30s " % items[item],end='')
  for t in experiment['trial']:
    print("& %10.5f " % t['msmout'][item], end='')
  print("\\\\")
    
  
