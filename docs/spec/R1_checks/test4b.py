import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import cmod
from test4 import prof, hs, shake_tf
from scipy.optimize import brentq
def column(hl,props,w,c):
    n=len(hl)+1; K=np.zeros((n,n),complex); Mlin=np.array([[2,1],[1,2]])/6.
    for j,(h,(Vs,rho,b,mm)) in enumerate(zip(hl,props)):
        Mm={'consistent':Mlin,'lumped':np.eye(2)/2,'mixed':0.5*Mlin+0.25*np.eye(2)}[mm]
        G=rho*Vs**2*cmod(b); K[j:j+2,j:j+2]+= G/h*np.array([[1,-1],[-1,1]]) - w**2*rho*h*Mm
    K[-1,-1]+=1j*w*c
    Fb=np.zeros(n,complex); Fb[-1]=1j*w*c*2.0
    return np.linalg.solve(K,Fb)
def tlm_outcrop(f,nhs=20,frac=0.1,mass='consistent'):
    w=2*np.pi*f; hl=[];props=[]
    for h,Vs,rho,b in prof:
        n=max(1,int(np.ceil(h/(frac*Vs/f)))); hl+=[h/n]*n; props+=[(Vs,rho,b,mass)]*n
    Vh,rh,bh=hs; D=1.5*Vh/f; h1=hl[-1]
    if nhs*h1<D:
        q=brentq(lambda q: h1*(q**nhs-1)/(q-1)-D,1.0000001,10); hh=[h1*q**i for i in range(nhs)]
    else: hh=[D/nhs]*nhs
    c=rh*Vh*np.sqrt(cmod(bh))
    U1=column(hl+hh,props+[(Vh,rh,bh,mass)]*nhs,w,c)
    U2=column(hh,[(Vh,rh,bh,mass)]*nhs,w,c)
    return U1[0]/U2[0]
if __name__ == '__main__':
    print('outcrop TF: TLM(soil+HS layers+dashpot) / TLM(HS layers+dashpot only) vs SHAKE surface/(2E_hs)')
    for f in [0.5,1.,1.7,2.5,4.,7.,10.,15.]:
        s=shake_tf(f); row=[]
        for nhs in [5,10,20]:
            t=tlm_outcrop(f,nhs=nhs,frac=0.05); row.append('n=%2d %.4f/%5.1fdeg'%(nhs,abs(t)/abs(s),np.degrees(np.angle(t/s))))
        tm=tlm_outcrop(f,nhs=20,frac=0.2,mass='mixed'); tc=tlm_outcrop(f,nhs=20,frac=0.2,mass='consistent')
        print('f=%5.2f |S|=%.3f  '%(f,abs(s))+'  '.join(row)+'  | h=lam/5 mixed %.3f consistent %.3f'%(abs(tm)/abs(s),abs(tc)/abs(s)))
