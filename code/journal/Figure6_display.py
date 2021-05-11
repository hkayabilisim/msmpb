from msmpb_wrapper import *
import pandas as pd
import json
import matplotlib.pyplot as plt

with open('../../results/Figure6.json', 'r') as inpfile:
    experiment = json.load(inpfile)

results = pd.DataFrame(columns=['L','nu','M','a0','cost','error'])
for t in experiment['trial']:
    row = pd.DataFrame({'L'       :[t['L']],
                        'nu'      :[t['nu']],
                        'a0'      :[t['a0']],
                        'M'       :[t['M']],
                        'cost'    :[t['msmout']['time_total']-t['msmout']['time_nlist']],
                        'error'   :[t['msmout']['deltaF/Fref']]
                       })
    results = results.append(row,ignore_index=True)
print(results)   
fig = plt.figure(figsize=(12, 4))
spec = fig.add_gridspec(ncols=2, nrows=1)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7)

for subplotindex in [0,1]:
    a0 = 9 if subplotindex == 0 else 12
    ax = fig.add_subplot(spec[0,subplotindex])
    
    idxL1 = (results['a0'] == a0) & (results['L'] == 1)
    idxL2 = (results['a0'] == a0) & (results['L'] == 2)
    idxL3 = (results['a0'] == a0) & (results['L'] == 3)
    
    ax.scatter(results[idxL1]['cost'], 1e4*results[idxL1]['error'], 
        label=r'L=1', edgecolor='k', marker='o')
    ax.scatter(results[idxL2]['cost'], 1e4*results[idxL2]['error'], 
        label=r'L=2', edgecolor='k', marker='s')
    ax.scatter(results[idxL3]['cost'], 1e4*results[idxL3]['error'], 
        label=r'L=3', edgecolor='k', marker='d')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("CPU time (sec)")
    ax.set_ylabel(r'relative force error x $10^4$')
    ax.set_title(r'$a_0={}$'.format(a0))
    ax.legend()
    ax.grid(b=True, which='both',linestyle='dotted',color='lightgrey')
    if subplotindex == 0:
      plt.yticks([1,2,3,4,6,10],['1','2','3','4','6','10'])    
      plt.xticks([0.3,0.4,0.6,1,2,3],['0.3','0.4','0.6','1','2','3'])
    else:
      plt.yticks([1,2,3,4,6,10],['1','2','3','4','6','10'])    
      plt.xticks([0.3,0.4,0.6,1,2,3],['0.3','0.4','0.6','1','2','3'])
    ax.legend(loc='lower left')
plt.savefig('../../results/Figure6.pdf')
