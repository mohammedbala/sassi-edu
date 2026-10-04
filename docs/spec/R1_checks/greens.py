import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import *
from scipy.special import hankel2, jv
def dH(n,x):  # d/dx H_n^(2)(x)
    return 0.5*(hankel2(n-1,x)-hankel2(n+1,x))
def green_point(md, m, n, rho_):
    """Kausel (1981) Eq.88-89, unit load at interface n, observed at interface m.
    returns dict: horizontal load (cos-theta comps) u_rho~, u_theta~ (with u_theta=-v~ sin), u_z~ ; vertical load u_rho, u_z"""
    kr,phx,phz,kl,phy=md
    xR=kr*rho_; xL=kl*rho_
    a=1/(4j)
    ur_h = a*( np.sum(phx[m]*phx[n]*kr*dH(1,xR)/kr) + (1/rho_)*np.sum(phy[m]*phy[n]*hankel2(1,xL)/kl) )
    ut_h = a*( (1/rho_)*np.sum(phx[m]*phx[n]*hankel2(1,xR)/kr) + np.sum(phy[m]*phy[n]*kl*dH(1,xL)/kl) )
    uz_h = -a*np.sum(phz[m]*phx[n]*hankel2(1,xR))
    ur_v = a*np.sum(phx[m]*phz[n]*hankel2(1,xR))
    uz_v = a*np.sum(phz[m]*phz[n]*hankel2(0,xR))
    return ur_h, ut_h, uz_h, ur_v, uz_v
