from msmpb_util import *
import pandas as pd

data = np.transpose(np.array([[  4, 4 ,  4,  4,  4,  4,  4,  4,  4,  4,  6,  6,  6,  6,  6,  6,  6,  6,  6], \
                 [ 108, 90, 80, 75, 72, 70, 64, 60, 50, 44, 90, 80, 75, 72, 70, 64, 60, 50, 44],\
                 [7.00,7.90,8.55,8.93,9.17,9.35,9.92,10.36,11.70,12.74,6.00,6.59,6.94,7.17,7.34,7.88,8.30,9.60,10.64]]))
box44 = {'name':'Box44'    ,'N':65850,'A':88}
results = np.zeros((data.shape[0],4))

fp = open("../results/Figure4-data.txt","w")
print("%2s %3s %5s %8s %8s %8s %8s"%("nu","M","a0","timeMSM","errMSM","timePME","errPME"),file=fp)
for i in range(data.shape[0]):
    nu = data[i][0]
    M = data[i][1]
    a0 = data[i][2]
    msmout = runMSM(box44,M,a0,nu)
    results[i][0]=msmout['time_total']
    results[i][1]=msmout['deltaF/Fref']
    
    pmeout = runPME(box44,M,a0,nu)
    results[i][2]=pmeout['time_total']
    results[i][3]=pmeout['deltaF/Fref']
    print("%2d %3d %5.2f %8.3f %8.3e %8.3f %8.3e"%(nu,M,a0,results[i][0],results[i][1],results[i][2],results[i][3]),file=fp)  
fp.close()


print("----- TABLE III, IV ------")
for order in [4,6]:
    print("Order is %d"%order)
    print("%20s %3s "%("","h"),end='')
    nuidx = data[:,0] == order
    for grid in np.unique(data[nuidx,1]):
        print("%5.2f "%(box44["A"]/grid),end='')
    print("")
    print("%20s %3s "%("Force Errorx1e4","PME"),end='')
    for grid in np.unique(data[nuidx,1]):
        idx = nuidx & (grid == data[:,1])    
        print("%5.2f "%(results[idx,3]*1e4),end='')
    print("")
    print("%20s %3s "%("Force Errorx1e4","MSM"),end='')
    for grid in np.unique(data[nuidx,1]):
        idx = nuidx & (grid == data[:,1])    
        print("%5.2f "%(results[idx,1]*1e4),end='')
    print("")
    print("%20s %3s "%("Time","PME"),end='')
    for grid in np.unique(data[nuidx,1]):
        idx = nuidx & (grid == data[:,1])    
        print("%5.2f "%(results[idx,2]),end='')
    print("")
    print("%20s %3s "%("Time","MSM"),end='')
    for grid in np.unique(data[nuidx,1]):
        idx = nuidx & (grid == data[:,1])    
        print("%5.2f "%(results[idx,0]),end='')
    print("")        

data = pd.read_csv('../results/Figure4-data.txt',delim_whitespace=True)

fig = plt.figure(figsize=(12, 4))
spec = fig.add_gridspec(ncols=2, nrows=1)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7)

theorymins = [1.1e-3,6e-4]
for i in [0,1]:
    order = 2*(i+2)
    idx = order == data['nu']

    xmin = 0.9*min(data[idx].min(axis=0)[[3,5]])
    xmax = 1.1*max(data[idx].max(axis=0)[[3,5]])
    ymin = 0.9*min(data[idx].min(axis=0)[[4,6]])
    ymax = 1.1*max(data[idx].max(axis=0)[[4,6]])

    ax = fig.add_subplot(spec[0,i])
    ax.scatter(data['timeMSM'][idx],data['errMSM'][idx], \
                        label='MSM', \
                        edgecolor='k')
    ax.scatter(data['timePME'][idx],data['errPME'][idx], \
                        label='PME', \
                        edgecolor='k',marker='s')
    theorymin =   theorymins[i]
    theorymax =   theorymin*(xmax/xmin)**(-order/3)
    ax.plot([xmin,xmax],[theorymin,theorymax])
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("Cost")
    ax.set_ylabel(r'Relative Force Error x $10^4$')
    ax.legend(loc='lower right')
    if order == 6:
        plt.xticks([1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8],['1.2','','','','2','','','','2.8'])
        plt.yticks([2e-4,3e-4,4e-4,6e-4,1e-3],['2','3','4','6','10'])
    if order == 4:
        plt.xticks([2,3,4],['2','3','4'])
        plt.yticks([4e-4,6e-4,1e-3],['4','6','10'])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.grid(b=True, which='both',linestyle='dotted',color='lightgrey')
    #ax.ticklabel_format(useOffset=False)
    #ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
        
    plt.title(r'$\nu='+str(order)+'$')
plt.savefig('../results/Figure4.pdf')


    
