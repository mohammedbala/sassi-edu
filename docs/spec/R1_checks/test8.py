import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import assemble, modes
from point3 import solve_point, Psi
layers = [(0.5, 150, 300, 1.9, 0.05, 0.05)]*8 + [(1.0, 250, 500, 2.0, 0.04, 0.04)]*10 + [(2.0, 400, 800, 2.1, 0.03, 0.03)]*10
om=2*np.pi*5.0; md=modes(assemble(layers),om); n=md[1].shape[0]; R0=0.9
sol={}
for ifc in [0,2,6]:
    for mu in [0,1]:
        sol[(ifc,mu)]=solve_point(layers,md,mu,R0,om,ifc,'x' if mu else 'z')[2]
def comps(ifc_load,mu,ifc_obs,r):
    P=Psi(md,mu,r)@sol[(ifc_load,mu)]
    return P[ifc_obs],P[n+ifc_obs],P[2*n+ifc_obs]
def block(i_ifc,j_ifc,dx,dy):
    """3x3 flexibility: displacement at node i (ifc i_ifc) due to unit loads at node j (ifc j_ifc); (dx,dy)=xi-xj"""
    r=np.hypot(dx,dy); c,s=dx/r,dy/r
    u,v,w=comps(j_ifc,1,i_ifc,r)      # horizontal load: u_rho=u cos, u_theta=-v sin, u_z=w cos
    uz_r,_,wz=comps(j_ifc,0,i_ifc,r)  # vertical load: u_rho=uz_r, u_z=wz
    return np.array([[u*c*c+v*s*s,(u-v)*s*c,uz_r*c],[(u-v)*s*c,u*s*s+v*c*c,uz_r*s],[w*c,w*s,wz]])
for (a,b,dx,dy) in [(0,6,2.7,1.8),(2,6,-3.1,4.4),(0,2,5.0,-1.0)]:
    Gij=block(a,b,dx,dy); Gji=block(b,a,-dx,-dy)
    print('pair ifc %d,%d: max|Gij-Gji^T|/max|Gij| = %.2e'%(a,b,np.abs(Gij-Gji.T).max()/np.abs(Gij).max()))
