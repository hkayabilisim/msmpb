from msmpb_wrapper import *  
        
data = [{'name':'spceN300' ,'N':300  ,'A':20},             
        {'name':'spceN600' ,'N':600  ,'A':20},             
        {'name':'spceN900' ,'N':900  ,'A':20},             
        {'name':'spceN2250','N':2250 ,'A':30}]             
        
results = []                                               
for benchmark in data:                                     
    for M in [24, 27, 30, 32, 35]:                           
        for a0 in [9,10,11]:                                   
            for nu in [4,6,8,10]:                                
                out = runMSM(benchmark,M,a0,nu)                    
                results.append(out)   
                


N = np.zeros((len(results),1))
a0 = np.zeros((len(results),1))
nu = np.zeros((len(results),1))
M = np.zeros((len(results),1))
error = np.zeros((len(results),1))
print("%4s %3s %8s %8s %2s %10s"%("N","M","a0","M/a0","nu","error"))
for i in range(len(results)):
    N[i] = results[i]["N"]
    a0[i] = results[i]["cutoff"]
    error[i] = results[i]["deltaF/Fref"]
    M[i] = results[i]["TopLevelMx"]
    nu[i] = results[i]["nu"]
    print("%4d %3d %8.3f %8.3f %2d %10.3e"%(N[i],M[i],a0[i],M[i]/a0[i],nu[i],error[i]))



fig, axes = plt.subplots(nrows=1, ncols=4)
plt.rcParams.update({'font.size': 12})
fig.set_size_inches(12, 9)
for order in [4,6,8,10]:
    ax = axes[order//2-2]
    for nparticles in np.unique(N):
        for cutoff in np.unique(a0):
            idx = (nu == order) & (N == nparticles) & (cutoff == a0);
            ax.loglog(M[idx]/(N[idx])**(1/3),error[idx],'gray',linewidth=2)
            
            
    nuidx = order == nu;
    MN = M[nuidx]/(N[nuidx]**(1/3))
    xmin=np.min(MN);
    xmax=np.max(MN);
    slope1minusNU = (10**(-order/2-1) )* ([xmin, xmax]/xmin)**(1-order);
    ax.loglog([xmin, xmax],slope1minusNU,'red',linewidth=2,label=r'$(M/N^{1/3})^{1-\nu}$');
    ax.set_ylim((1e-9,1e-3))  
    ax.set_xticks([2, 3,4])
    ax.set_xticklabels(['2','3','4'])
    ax.set_ylabel('relative force error')
    ax.set_xlabel(r'$M/N^{1/3}$')
    ax.legend(loc='upper right')
    ax.set_title(r'$\nu=%d$'%order)
    ax.grid(b=True, which='both',linestyle='dotted',color='lightgrey')
fig.tight_layout(pad=1.0)       
fig.savefig('../../results/Figure1.pdf', dpi=100)

