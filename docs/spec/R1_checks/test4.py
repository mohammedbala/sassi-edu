import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import cmod
# SHAKE exact (Schnabel et al 1972 / SHAKE91 CXSOIL+AMP): surface / outcrop-at-halfspace
prof=[(10.,150.,1.8,0.05),(15.,250.,1.9,0.04),(20.,400.,2.0,0.03)]  # h,Vs,rho,beta
hs=(800.,2.2,0.01)
def shake_tf(f):
    w=2*np.pi*f
    E=1.;F=1.
    lay=prof+[(0.,)+hs]
    for i in range(len(prof)):
        h,Vs,rho,b=prof[i]; G=rho*Vs**2*cmod(b); V=np.sqrt(G/rho)
        h2,Vs2,rho2,b2=lay[i+1]; G2=rho2*Vs2**2*cmod(b2)
        a=np.sqrt(rho*G/(rho2*G2))
        ex=np.exp(1j*w*h/V)
        E,F=0.5*(1+a)*E*ex+0.5*(1-a)*F/ex, 0.5*(1-a)*E*ex+0.5*(1+a)*F/ex
    return 2.0/(2*E)   # surface (E1+F1=2) over outcrop 2E_hs
def tlm_tf(f,nhs=20,hmax_frac=0.1,mass='consistent'):
    w=2*np.pi*f
    hl=[];props=[]
    for h,Vs,rho,b in prof:
        n=max(1,int(np.ceil(h/(hmax_frac*Vs/f)))); hl+= [h/n]*n; props+=[(Vs,rho,b)]*n
    Vh,rh,bh=hs; lam=Vh/f; D=1.5*lam
    # variable-depth halfspace: geometric progression starting at last soil sublayer thickness
    h1=hl[-1]; 
    from scipy.optimize import brentq
    if nhs*h1<D:
        q=brentq(lambda q: h1*(q**nhs-1)/(q-1)-D,1.0000001,10)
        hh=[h1*q**i for i in range(nhs)]
    else: hh=[D/nhs]*nhs
    hl+=hh; props+=[(Vh,rh,bh)]*nhs
    n=len(hl)+1
    K=np.zeros((n,n),complex)
    Mlin=np.array([[2,1],[1,2]])/6.
    Mm={'consistent':Mlin,'lumped':np.eye(2)/2,'mixed':0.5*Mlin+0.25*np.eye(2)}[mass]
    for j,(h,(Vs,rho,b)) in enumerate(zip(hl,props)):
        G=rho*Vs**2*cmod(b)
        K[j:j+2,j:j+2]+= G/h*np.array([[1,-1],[-1,1]]) - w**2*rho*h*Mm
    Vstar=Vh*np.sqrt(cmod(bh)); c=rh*Vstar
    K[-1,-1]+=1j*w*c
    Eb=1.0; Fb=np.zeros(n,complex); Fb[-1]=1j*w*c*2*Eb
    U=np.linalg.solve(K,Fb)
    kstar=w/Vstar; Dtot=sum(hh)
    Etop=Eb*np.exp(-1j*kstar*Dtot)
    return U[0]/(2*Etop)
def shake_within(f):
    w=2*np.pi*f; E=1.;F=1.; lay=prof+[(0.,)+hs]
    for i in range(len(prof)):
        h,Vs,rho,b=prof[i]; G=rho*Vs**2*cmod(b); V=np.sqrt(G/rho)
        h2,Vs2,rho2,b2=lay[i+1]; G2=rho2*Vs2**2*cmod(b2); a=np.sqrt(rho*G/(rho2*G2)); ex=np.exp(1j*w*h/V)
        E,F=0.5*(1+a)*E*ex+0.5*(1-a)*F/ex, 0.5*(1-a)*E*ex+0.5*(1+a)*F/ex
    return 2.0/(E+F)
def tlm_within(f,nhs=20,frac=0.1,mass='consistent',uniform=False):
    w=2*np.pi*f
    hl=[];props=[]
    for h,Vs,rho,b in prof:
        n=max(1,int(np.ceil(h/(frac*Vs/f)))); hl+= [h/n]*n; props+=[(Vs,rho,b)]*n
    nsoil=len(hl)
    Vh,rh,bh=hs; D=1.5*Vh/f; h1=hl[-1]
    from scipy.optimize import brentq
    if (not uniform) and nhs*h1<D:
        q=brentq(lambda q: h1*(q**nhs-1)/(q-1)-D,1.0000001,10); hh=[h1*q**i for i in range(nhs)]
    else: hh=[D/nhs]*nhs
    hl+=hh; props+=[(Vh,rh,bh)]*nhs
    n=len(hl)+1; K=np.zeros((n,n),complex)
    Mlin=np.array([[2,1],[1,2]])/6.; Mm={'consistent':Mlin,'lumped':np.eye(2)/2,'mixed':0.5*Mlin+0.25*np.eye(2)}[mass]
    for j,(h,(Vs,rho,b)) in enumerate(zip(hl,props)):
        G=rho*Vs**2*cmod(b); K[j:j+2,j:j+2]+= G/h*np.array([[1,-1],[-1,1]]) - w**2*rho*h*Mm
    c=rh*Vh*np.sqrt(cmod(bh)); K[-1,-1]+=1j*w*c
    Fb=np.zeros(n,complex); Fb[-1]=1j*w*c*2.0
    U=np.linalg.solve(K,Fb)
    return U[0]/U[nsoil], hh
if __name__ == '__main__':
    print('--- within-motion TF: surface / halfspace-top (independent of extra-layer phase)')
    for f in [0.5,1.,1.7,2.5,4.,7.,10.,15.]:
        s=shake_within(f)
        out=[]
        for nhs in [5,10,20]:
            t,hh=tlm_within(f,nhs=nhs)
            out.append('n=%d: %.4f/%.1fdeg'%(nhs,abs(t)/abs(s),np.degrees(np.angle(t/s))))
        t2,_=tlm_within(f,nhs=20,frac=0.05)
        print('f=%5.2f |SHAKE|=%.4f  '%(f,abs(s))+'  '.join(out)+'   fine soil n=20: %.4f/%.1fdeg'%(abs(t2)/abs(s),np.degrees(np.angle(t2/s))))
