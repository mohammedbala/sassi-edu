import numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.special import hankel2, jv
from scipy.optimize import brentq
from tlm import assemble, modes
from point3 import solve_point
layers = [(0.5, 150, 300, 1.9, 0.05, 0.05)]*8 + [(1.0, 250, 500, 2.0, 0.04, 0.04)]*10 + [(2.0, 400, 800, 2.1, 0.03, 0.03)]*10
f=5.0; om=2*np.pi*f
md=modes(assemble(layers),om); kr,phx,phz,kl,phy=md
def I3_at_R(k,a): return np.pi/(2j*k**2)*jv(1,k*a)*hankel2(1,k*a)-1/(2*k**2)
def disk_avg_x(m,n,a):
    q=1/(np.pi*a**2)
    return q*(np.sum(phx[m]*phx[n]*I3_at_R(kr,a))+np.sum(phy[m]*phy[n]*I3_at_R(kl,a)))
def disk_avg_z(m,n,a):
    q=1/(np.pi*a**2)
    return 2*q*np.sum(phz[m]*phz[n]*I3_at_R(kr,a))
for R0 in [0.5,0.9,1.5]:
    Ux,_,_=solve_point(layers,md,1,R0,om,0,'x'); Uz,_,_=solve_point(layers,md,0,R0,om,0,'z')
    ax=brentq(lambda a: abs(disk_avg_x(0,0,a))-abs(Ux[0]),0.05,5)
    az=brentq(lambda a: abs(disk_avg_z(0,0,a))-abs(Uz[0]),0.05,5)
    print('R0=%.2f  FE-core axis ux=%s uz=%s | equivalent uniform-disk radius (avg displ.) a_x/R0=%.3f a_z/R0=%.3f'%(R0,np.round(Ux[0],9),np.round(Uz[0],9),ax/R0,az/R0))
print('equal-area radius for square mesh h: a=h/sqrt(pi)=0.564h ; with R0=0.9h -> a/R0=%.3f'%(1/np.sqrt(np.pi)/0.9))
