import numpy as np
from tlm import cmod
# exact 2DOF hysteretic base-excited system: total-acceleration TF of DOF 1 (top)
m1,m2=1.0,1.5; k1,k2=(2*np.pi*3.)**2, (2*np.pi*5.)**2*2
b=0.05
def tf(w):
    K=np.array([[k1,-k1],[-k1,k1+k2]])*cmod(b); M=np.diag([m1,m2])
    r=np.array([1.,1.])
    # relative u: (K - w^2 M) u = w^2 M r ag (ag displacement=1 -> abs = u + r)
    u=np.linalg.solve(K-w**2*M, w**2*M@r)
    return u+r
def fit5(ws,Us):
    A=np.array([[w**4,w**2,1,-w**2*U,-U] for w,U in zip(ws,Us)])
    rhs=np.array([w**4*U for w,U in zip(ws,Us)])
    return np.linalg.solve(A,rhs)
def evalc(c,w): return (c[0]*w**4+c[1]*w**2+c[2])/(w**4+c[3]*w**2+c[4])
ws=2*np.pi*np.array([1.0,2.5,4.0,6.0,9.0])
for dof in [0,1]:
    Us=np.array([tf(w)[dof] for w in ws]); c=fit5(ws,Us)
    wt=2*np.pi*np.linspace(0.5,12,300)
    err=max(abs(evalc(c,w)-tf(w)[dof])/abs(tf(w)[dof]) for w in wt)
    print('dof',dof,'max rel error of 5-point 2DOF interpolant over 0.5-12 Hz:',err)
# scaled version (better conditioning): w -> w/wref
