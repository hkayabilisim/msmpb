from msmpb_util import *
import pandas as pd
import json
import matplotlib.pyplot as plt

with open('../../results/Figure5.json', 'r') as inpfile:
    experiment = json.load(inpfile)

results = pd.DataFrame(columns=['nu','M','a0','msm_cost','pme_cost','msm_uerr','pme_uerr'])
all_time = 0
for t in experiment['trial']:
    msm_cost = 0
    msm_uerr = 0
    for snapshot in t['msmout']:
        msm_cost = msm_cost + snapshot['time_total']-snapshot['time_nlist']
        msm_uerr = msm_uerr + snapshot['potabserror']
        all_time += snapshot['time_total']
    msm_cost = msm_cost / len(t['msmout'])
    msm_uerr = msm_uerr / len(t['msmout'])
    
    pme_cost = 0
    pme_uerr = 0
    for snapshot in t['pmeout']:
        pme_cost = pme_cost + snapshot['time_total']-snapshot['time_nlist']
        pme_uerr = pme_uerr + snapshot['potabserror']
        all_time += snapshot['time_total']
    pme_cost = pme_cost / len(t['pmeout'])
    pme_uerr = pme_uerr / len(t['pmeout'])
    
    row = pd.DataFrame({'nu'      :[t['nu']],
                        'a0'      :[t['a0']],
                        'M'       :[t['M']],
                        'msm_cost':[msm_cost],
                        'msm_uerr':[msm_uerr],                        
                        'pme_cost':[pme_cost],
                        'pme_uerr':[pme_uerr],
                       })
    results = results.append(row,ignore_index=True)

print(results)
print(all_time)
factor = 332.0636;
fig = plt.figure(figsize=(12, 4))
spec = fig.add_gridspec(ncols=2, nrows=1)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7)

for subplotindex in [0,1]:
    a0 = 9 if subplotindex == 0 else 12
    ax = fig.add_subplot(spec[0,subplotindex])
    
    idx = results['a0'] == a0

    ax.scatter(results[idx]['msm_cost'], factor*results[idx]['msm_uerr'], 
        label=r'MSM', edgecolor='k', marker='o')
    ax.scatter(results[idx]['pme_cost'], factor*results[idx]['pme_uerr'], 
        label=r'PME', edgecolor='k', marker='s')
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("CPU time (sec)")
    ax.set_ylabel(r'avg. abs. potential error (kcal/mol)')
    ax.set_title(r'$a_0={}$'.format(a0))
    ax.legend()
    ax.grid(b=True, which='both',linestyle='dotted',color='lightgrey')
    plt.yticks([1e-2,1e-1,1e0,1e1,1e2],['0.01','0.1','1','10','100'])    
    plt.xticks([0.3,0.4,0.6,1],['0.3','0.4','0.6','1'])
    if subplotindex == 1:
        plt.xticks([0.6,0.8,1.2,2],['0.6','0.8','1.2','2'])
        ax.set_xlim(.6,2)
    ax.legend(loc='lower left')

plt.savefig('../../results/Figure5.pdf')
