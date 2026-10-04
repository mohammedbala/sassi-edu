import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import assemble, modes, direct_flex
from scipy.integrate import quad
layers = [(1.0, 150, 300, 1.9, 0.05, 0.05)]*6 + [(2.0, 300, 600, 2.0, 0.05, 0.05)]*6
f=4.0; om=2*np.pi*f
mats=assemble(layers); md=modes(mats,om); kr,phx,phz,kl,phy=md; n=phx.shape[0]
def eq47(m,nn,x):
    E=np.exp(-1j*kr*abs(x)); s=1 if x>=0 else -1
    uxx=np.sum(phx[m]*phx[nn]*E/kr)/(2j)
    uxz=s*np.sum(phx[m]*phz[nn]*E/kr)/(2j)
    uzx=s*np.sum(phz[m]*kr**0*phx[nn]*E/kr)/(2j)
    uzz=np.sum(phz[m]*phz[nn]*E/kr)/(2j)
    return uxx,uxz,uzx,uzz
def numeric(m,nn,x,kmax=60.):
    def F(k): return direct_flex(mats,om,k)
    out=[]
    for (a,b) in [(m,nn),(m,n+nn),(n+m,nn),(n+m,n+nn)]:
        re=quad(lambda k: (F(k)[a,b]*np.exp(-1j*k*x)).real,-kmax,kmax,limit=2000,points=[0])[0]
        im=quad(lambda k: (F(k)[a,b]*np.exp(-1j*k*x)).imag,-kmax,kmax,limit=2000,points=[0])[0]
        out.append((re+1j*im)/(2*np.pi))
    return out
for x in [3.0,-3.0,8.0]:
    a=eq47(0,3,x); b=numeric(0,3,x)
    print('x=%5.1f'%x,' Eq47:',np.round(a,8)); print('        numeric:',np.round(b,8))
