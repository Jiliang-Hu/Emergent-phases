import numpy as np, pandas as pd, sys, matplotlib.figure as mpfig, matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


EPS=10.**-8 ### Numerical threshold
dispersal=10**-6

def run(
        S=20, #Species number
        cp=100., #Coupling strength
        days=15,
        pH0=7,
        dilution=0.1, #if >0, daily dilution
        pc=2.5, #gate width
        kgrowth=10., #growth rate
        **kwargs
            ):

    ppref=kwargs.get('ppref',np.random.uniform(pH0-pc,pH0+pc,S))
    cs= kwargs.get('cs',np.random.uniform(-1,1,S)*cp )
    x0=kwargs.get('x0',np.concatenate([[pH0],np.log(np.random.random(S)*0.01)]) )

    def fun(p):
        ingate=(np.abs(p-ppref)<pc ).astype('float')
        return  ingate 

    def eqs(t,x):
        p=x[0]
        x=np.exp(np.clip(x[1:],np.log(EPS),0))
        dx= x * (1-x)* kgrowth*fun(p) 
        if dilution<0:
            #Continuous-time dilution
            dx+= np.log(-dilution)*x                
        dx+= dispersal
        dp=np.sum(cs*x) 
        if dilution<0:
            dp+= np.log(-dilution)*(p-pH0)
        if p>14:
            dp=min(0,dp)
        if p<0:
            dp=max(0,dp)
        return np.concatenate([[dp],dx/x])
    
    time=[]
    pH=[]
    Ns=[]
    for day in range(days):
        sol=solve_ivp(fun=eqs, t_span=(0,1.),y0=x0)
        pH=np.concatenate([pH,sol.y[0]])
        if len(Ns):
            Ns=np.concatenate([Ns,np.exp(sol.y[1:].T)])
        else:
            Ns=np.exp(sol.y[1:].T)
        time=np.concatenate([time,sol.t+day])
        x0=sol.y[:,-1].copy()
        if dilution>0:
            #Discrete daily dilution
            x0[0]=(1-dilution)*pH0 + dilution* x0[0]
            x0[1:]=np.clip(x0[1:]+np.log(dilution),np.log(EPS),None)
    return time[:-1],pH[:-1],Ns[:-1]


def measure(ts,pH,Ns):
    Nf=Ns[-1]
    S=len(Nf)
    days=sorted(set(ts//1))
    
    observation_period=int(max(days))//3
    Nmean=np.array([ np.mean(Ns[ts//1 ==z],axis=0) for z in days ])
    nrel=Nmean/np.sum(Nmean,axis=1).reshape((-1,1))
    pHmean=np.array([ np.mean(pH[ts//1 ==z]) for z in days ])
    Ntot=np.sum(Nmean,axis=1)

    alive=np.where(Nmean[-1]>dispersal*S)[0]
    nalive=len(alive)

    CV=np.std(Nmean[-observation_period:,alive],axis=0)/(0.01+np.mean(Nmean[-observation_period:,alive],axis=0))
    meanCV=np.mean(CV)
    fluct=(meanCV > 0.25)
    
    return {'alive':float(nalive)/S,'fluct':float(fluct),}





Ss=[1,2,3,4,6,8,10,12,15,18,21,24,30,36,42,48]
cps=np.logspace(0,4,13)
replicas=500
days=20
table=[]
dilution=-0.033

fname='pHdata.csv'
df=None
for S in Ss:
    for cp in cps:
        print 'RUN:', S, cp
        nrep=replicas
        if S==Ss[-1]:
            nrep=replicas*2
        for replica in range(nrep):
            print '         ',replica
            opts={'days':days,'cp':cp,'dilution':dilution,}
            ts,pH,Ns=run(S=S,**opts)  
            dic={'S':S, 'cp':cp,'replica':replica}
            dic.update(measure(ts,pH,Ns))
            table.append(dic)
df=pd.DataFrame(table)  
df.to_csv(fname)
