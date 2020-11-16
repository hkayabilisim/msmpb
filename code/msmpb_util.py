import subprocess
import numpy as np                                         
import matplotlib.pyplot as plt                                
import subprocess
import json

def calculateM(eta,N):
    return int(np.ceil(pow(eta * N, 1/3)))

def getFactors(M):
    factors = []
    divisor = 2
    while M != 1:
        if M % divisor == 0:
            M = M/divisor
            if divisor not in factors:
                factors.append(divisor)
        else:
            divisor = divisor + 1
    return factors

def isFFTFriendly(M):
    factors = getFactors(M)
    if (2 not in factors) and (3 not in factors) and (5 not in factors) and (7 not in factors):
        return False
    else:
        return True

def makeFFTFriendly(M,decrease=True):
    if decrease == True:
        step = -1
    else:
        step =  1
    while not isFFTFriendly(M):
        M = M + step
    return M

def nextFFTFriendly(M):  
    while not isFFTFriendly(M+1):
        M = M + 1
    return M + 1

def prevFFTFriendly(M):  
    while not isFFTFriendly(M-1):
        M = M - 1
    return M - 1

def runMSM(benchmark,M,a0,nu,L=1):
    A = benchmark['A']
    name = benchmark['name'];
    N = benchmark['N']
    h = A / M;
    nbar = 2 * a0 / h;
    out = subprocess.run(['./msmpb','../data/'+name,'-L',str(L),'--a0',str(a0),'--nu',str(nu),'-M',str(M)]\
                     ,stdout=subprocess.PIPE)
    return json.loads(out.stdout.decode("utf-8"))

def runPME(benchmark,M,a0,nu):
    name = benchmark['name'];
    out = subprocess.run(['./pme','../data/'+name,'--a0',str(a0),'--nu',str(nu),'-M',str(M)]\
                     ,stdout=subprocess.PIPE)
    return json.loads(out.stdout.decode("utf-8"))


def fudgeFactor(results):
    k = len(results)
    deltaF = np.zeros((k,1))
    deltaFest = np.zeros((k,1))
    nu = np.zeros((k,1))
    for i in range(0,k):
        deltaF[i] = results[i]["deltaF/Fref"]
        deltaFest[i] = results[i]["deltaFest/Fref"]
        nu[i] = results[i]["nu"]
    factors = []
    for order in [4,6,8,10]:
        nuidx = nu == order
        true = deltaF[nuidx]
        appr = deltaFest[nuidx]
        factor =  10**((np.sum(np.log10(true))-np.sum(np.log10(appr)))/len(true));
        factors.append(factor)
    return factors

def displayResult(results,fudgefactors):
    print("%20s %5s %2s %5s %2s %8s %2s %10s %10s %s"%("data","N","A","nbar","M","a0","nu","true","est","|log(true)-log(est)|"))
    for r in results:
        M = r['TopLevelMx']
        A = r['Edge row1'][0]
        nbar = 2 * r['cutoff'] * M / A
        nu = r['nu']
        est = r["deltaFest/Fref"] * fudgefactors[int(nu/2-2)]
        true = r["deltaF/Fref"]
        diff = np.abs(np.log(est)-np.log(true))
        print("%20s %5d %2d %5.2f %2d %8.3f %2d %10.3e %10.3e %5.3f" \
              %(r['data'],r['N'],A,nbar,M,r['cutoff'],nu,true,est,diff))


def plotEstimatedVersusTrue(ax,results,fudgefactors):
    k = len(results)
    estimated = np.zeros((k,1))
    true = np.zeros((k,1))
    nu = np.zeros((k,1))
    for i in range(k):
        estimated[i] = results[i]["deltaFest/Fref"]
        true[i] = results[i]["deltaF/Fref"]       
        nu[i] = results[i]["nu"]
    for i in [0,1,2,3]:
        order = 2*(i+2)
        idx = nu == order
        print("Fudge factor: %f"%(fudgefactors[i]))
        #*estimated[idx]
        ax.scatter(estimated[idx]*fudgefactors[i],true[idx],500, \
                    label=r'$\nu='+str(order)+'$', \
                    edgecolor='k')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("estimated")
    plt.ylabel("true")
    plt.xlim(1e-8, 1e-2)
    plt.ylim(1e-8, 1e-2)
    ax.plot([1e-8, 1e-2],[1e-8, 1e-2],label=None)
    #plt.grid(b=True,which="major",linestyle='dotted')
    ax.grid(linestyle='dotted')
    ax.legend()
    #plt.show()
    
def plotEstimatedVersusTrueWRTa0(ax,results,fudgefactors):
    k = len(results)
    estimated = np.zeros((k,1))
    true = np.zeros((k,1))
    a0 = np.zeros((k,1))
    for i in range(k):
        nu = results[i]["nu"]
        estimated[i] = fudgefactors[int(nu/2-2)]*results[i]["deltaFest/Fref"]
        true[i] = results[i]["deltaF/Fref"]       
        a0[i] = results[i]["cutoff"]
    for cutoff in np.unique(a0):
        idx = a0 == cutoff
        ax.scatter(estimated[idx],true[idx],500, \
                    label=r'$a_0='+str(cutoff)+'$', \
                    edgecolor='k')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("estimated")
    plt.ylabel("true")
    plt.xlim(1e-8, 1e-2)
    plt.ylim(1e-8, 1e-2)
    ax.plot([1e-8, 1e-2],[1e-8, 1e-2],label=None)
    #plt.grid(b=True,which="major",linestyle='dotted')
    ax.grid(linestyle='dotted')
    ax.legend()
    #plt.show()    
    
def plotEstimatedVersusTrueWRTM(ax,results,fudgefactors):
    k = len(results)
    estimated = np.zeros((k,1))
    true = np.zeros((k,1))
    M = np.zeros((k,1))
    for i in range(k):
        nu = results[i]["nu"]
        estimated[i] = fudgefactors[int(nu/2-2)]*results[i]["deltaFest/Fref"]
        true[i] = results[i]["deltaF/Fref"]       
        M[i] = results[i]["TopLevelMx"]
    for gs in np.unique(M):
        idx = M == gs
        ax.scatter(estimated[idx],true[idx],500, \
                    label=r'$M='+str(gs)+'$', \
                    edgecolor='k')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("estimated")
    plt.ylabel("true")
    plt.xlim(1e-8, 1e-2)
    plt.ylim(1e-8, 1e-2)
    ax.plot([1e-8, 1e-2],[1e-8, 1e-2],label=None)
    ax.grid(linestyle='dotted')
    ax.legend()
    
def plotEstimatedVersusTrueWRTN(ax,results,fudgefactors):
    k = len(results)
    estimated = np.zeros((k,1))
    true = np.zeros((k,1))
    N = np.zeros((k,1))
    for i in range(k):
        nu = results[i]["nu"]
        estimated[i] = fudgefactors[int(nu/2-2)]*results[i]["deltaFest/Fref"]
        true[i] = results[i]["deltaF/Fref"]       
        N[i] = results[i]["N"]
    for nn in np.unique(N):
        idx = N == nn
        ax.scatter(estimated[idx],true[idx],500, \
                    label=r'$N='+str(nn)+'$', \
                    edgecolor='k')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("estimated")
    plt.ylabel("true")
    plt.xlim(1e-8, 1e-2)
    plt.ylim(1e-8, 1e-2)
    ax.plot([1e-8, 1e-2],[1e-8, 1e-2],label=None)
    ax.grid(linestyle='dotted')
    ax.legend()

