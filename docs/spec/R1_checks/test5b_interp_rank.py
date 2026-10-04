"""5-point 2DOF rational interpolation on SDOF data: hysteretic SDOF -> rank-deficient (rank 4) but the
minimum-norm (lstsq) solution still interpolates exactly; viscous SDOF -> small non-zero error."""
import numpy as np
from tlm import cmod
def fit5(ws, Us, scale):
    x = ws/scale
    A = np.array([[w**4, w**2, 1, -w**2*U, -U] for w, U in zip(x, Us)]); rhs = np.array([w**4*U for w, U in zip(x, Us)])
    c, res, rk, sv = np.linalg.lstsq(A, rhs, rcond=1e-10)
    return c, rk, sv[0]/sv[-1]
def ev(c, w, scale): x = w/scale; return (c[0]*x**4 + c[1]*x**2 + c[2])/(x**4 + c[3]*x**2 + c[4])
w0 = 2*np.pi*4; ks = w0**2*cmod(0.05)
sd = lambda w: ks/(ks - w**2)
vis = lambda w: (w0**2 + 2j*0.05*w0*w)/(w0**2 - w**2 + 2j*0.05*w0*w)
ws = 2*np.pi*np.array([2., 3., 4., 5., 6.]); wt = 2*np.pi*np.linspace(2, 6, 400)
for name, fn in [('hysteretic SDOF', sd), ('viscous SDOF', vis)]:
    c, rk, cond = fit5(ws, fn(ws), ws[-1])
    err = max(abs(ev(c, w, ws[-1]) - fn(w))/abs(fn(w)) for w in wt)
    print('%-16s rank=%d cond=%.1e  max rel interp error in window = %.2e' % (name, rk, cond, err))
