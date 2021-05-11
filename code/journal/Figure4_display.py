from msmpb_wrapper import *
import pandas as pd
import json
import matplotlib.pyplot as plt

with open('../../results/Figure4.json', 'r') as inpfile:
    experiment = json.load(inpfile)

results = pd.DataFrame(columns=['nu','M','a0','msm_cost','pme_cost'])
for t in experiment['trial']:
    row = pd.DataFrame({'nu'      :[t['nu']],
                        'a0'      :[t['a0']],
                        'M'       :[t['M']],
                        'msm_cost':[t['msmout']['time_total']-t['msmout']['time_nlist']],
                        'msm_ferr':[t['msmout']['deltaF/Fref']],
                        'pme_cost':[t['pmeout']['time_total']-t['pmeout']['time_nlist']],
                        'pme_beta':[t['pmeout']['beta']],
                        'pme_ferr':[t['pmeout']['deltaF/Fref']],
                       })
    results = results.append(row,ignore_index=True)
print(results)   
fig = plt.figure(figsize=(12, 4))
spec = fig.add_gridspec(ncols=2, nrows=1)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7)

for subplotindex in [0,1]:
    a0 = 9 if subplotindex == 0 else 12
    ax = fig.add_subplot(spec[0,subplotindex])
    
    idx = results['a0'] == a0

    ax.scatter(results[idx]['msm_cost'], 1e4*results[idx]['msm_ferr'], 
        label=r'MSM', edgecolor='k', marker='o')
    ax.scatter(results[idx]['pme_cost'], 1e4*results[idx]['pme_ferr'], 
        label=r'PME', edgecolor='k', marker='s')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("CPU time (sec)")
    ax.set_ylabel(r'relative force error x $10^4$')
    ax.set_title(r'$a_0={}$'.format(a0))
    ax.legend()
    ax.grid(b=True, which='both',linestyle='dotted',color='lightgrey')
    plt.yticks([1,2,5,10,20],['1','2','5','10','20'])    
    plt.xticks([0.3,0.4,0.6,1],['0.3','0.4','0.6','1'])
    if subplotindex == 1:
        plt.xticks([0.6,0.8,1.2,2],['0.6','0.8','1.2','2'])
        ax.set_xlim(.6,2)
    ax.legend(loc='lower left')

plt.savefig('../../results/Figure4.pdf')
