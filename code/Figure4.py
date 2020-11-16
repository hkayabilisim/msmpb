from msmpb_util import *

data = np.transpose(np.array([[  4, 4 ,  4,  4,  4,  4,  4,  4,  4,  4,  6,  6,  6,  6,  6,  6,  6,  6,  6], \
                 [ 108, 90, 80, 75, 72, 70, 64, 60, 50, 44, 90, 80, 75, 72, 70, 64, 60, 50, 44],\
                 [7.00,7.90,8.55,8.93,9.17,9.35,9.92,10.36,11.70,12.74,6.00,6.59,6.94,7.17,7.34,7.88,8.30,9.60,10.64]]))
box44 = {'name':'Box44'    ,'N':65850,'A':88}
results = np.zeros((data.shape[0],4))

print("%2s %3s %5s %8s %8s %8s %8s"%("nu","M","a0","timeMSM","errMSM","timePME","errPME"))
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
    print("%2d %3d %5.2f %8.3f %8.3e %8.3f %8.3e"%(nu,M,a0,results[i][0],results[i][1],results[i][2],results[i][3]))

print("----- TABLE 4 ------")
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
    
plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 22})
for i in [1,2]:
    order = 2*(i+1)
    idx = order == data[:,0]

    xmin = 0.9*min(results[idx,:].min(axis=0)[[0,2]])
    xmax = 1.1*max(results[idx,:].max(axis=0)[[0,2]])
    ax = plt.subplot(1,2,i)
    ax.scatter(results[idx,0],results[idx,1],500, \
                        label='MSM', \
                        edgecolor='k')
    ax.scatter(results[idx,2],results[idx,3],500, \
                        label='PME', \
                        edgecolor='k',marker='s')
    theorymin =   1e-3
    theorymax =   results[0,1]*(xmax/xmin)**(-order/3)
    ax.plot([xmin,xmax],[theorymin,theorymax])
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Total Time")
    plt.ylabel("Error")
    plt.legend()
    plt.xlim([xmin,xmax])
    plt.ylim([0.9*min(results[idx,:].min(axis=0)[[1,3]]),1.1*max(results[idx,:].max(axis=0)[[1,3]])])
    plt.title(r'$\nu='+str(order)+'$')
plt.savefig('../results/Figure4.pdf')      
