import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import *
from scipy.special import hankel2
# ---- test 2: fundamental Rayleigh velocity, uniform deep stratum, nu=0.25
Vs=100.; nu=0.25; Vp=Vs*np.sqrt(2*(1-nu)/(1-2*nu)); rho=1.0
f=20.0; lam_s=Vs/f
for mass in ['consistent','lumped','mixed']:
    h=lam_s/20; nl=int(4*lam_s/h)
    layers=[(h,Vs,Vp,rho,1e-6,1e-6)]*nl
    mats=assemble(layers,mass=mass)
    om=2*np.pi*f
    kr,phx,phz,kl,phy=modes(mats,om)
    c=om/np.real(kr); cR=np.max(c[np.abs(np.imag(kr))<1e-6*np.abs(kr)]) if False else None
    # fundamental Rayleigh = largest real wavenumber among propagating modes
    prop=np.abs(np.imag(kr))<1e-3*np.abs(np.real(kr))
    kmax=np.max(np.real(kr[prop]))
    print(mass,'h=lambda/20: c_R/Vs =',om/kmax/Vs,' (exact 0.919402)')
    klmax=np.max(np.real(kl[np.abs(np.imag(kl))<1e-3*np.abs(np.real(kl))]))
    print('   Love max k -> c/Vs =',om/klmax/Vs,' (uniform stratum: -> 1.0)')
