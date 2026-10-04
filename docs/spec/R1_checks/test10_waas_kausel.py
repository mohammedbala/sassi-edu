"""Waas/SASSI layer form (z down, (A k^2 + i k B_W + G - w^2 M)) vs Kausel Table 1 form (z up, i*u_z scaled):
T^-1 L_W T == L_K with T = diag(1, i, 1, i), dof order (u1, v1, u2, v2) / (x1, z1, x2, z2)."""
import numpy as np
lam, G, rho, h, k, w = 3.0, 1.3, 1.7, 0.8, 0.7, 2.1
lp = lam + 2*G
A = h/6*np.array([[2*lp, 0, lp, 0], [0, 2*G, 0, G], [lp, 0, 2*lp, 0], [0, G, 0, 2*G]])
BW = 0.5*np.array([[0, -(lam-G), 0, (lam+G)], [(lam-G), 0, (lam+G), 0], [0, -(lam+G), 0, (lam-G)], [-(lam+G), 0, -(lam-G), 0]])
Gm = 1/h*np.array([[G, 0, -G, 0], [0, lp, 0, -lp], [-G, 0, G, 0], [0, -lp, 0, lp]])
M = rho*h/6*np.array([[2, 0, 1, 0], [0, 2, 0, 1], [1, 0, 2, 0], [0, 1, 0, 2]])
LW = A*k**2 + 1j*k*BW + Gm - w**2*M
BK = 0.5*np.array([[0, lam-G, 0, -(lam+G)], [lam-G, 0, lam+G, 0], [0, lam+G, 0, -(lam-G)], [-(lam+G), 0, -(lam-G), 0]])
LK = A*k**2 + k*BK + Gm - w**2*M
T = np.diag([1, 1j, 1, 1j])
print('max|T^-1 L_W T - L_K| =', np.abs(np.linalg.inv(T)@LW@T - LK).max())
print('B_W skew-symmetric:', np.allclose(BW, -BW.T), '  B_K symmetric:', np.allclose(BK, BK.T))
