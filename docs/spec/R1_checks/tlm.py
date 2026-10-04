"""Scratch prototype to verify thin-layer-method (TLM) formulas for the R1 theory note.
Kausel convention: exp(i w t - i k x), z UP, vector {u_x, u_y, i*u_z}; interfaces numbered from top.
"""
import numpy as np
np.seterr(all="ignore")
from scipy.linalg import eig, solve
from scipy.special import hankel2, jv

def cmod(beta):
    """SASSI/SHAKE/LUSH complex modulus factor (1-2b^2+2ib sqrt(1-b^2))."""
    return 1 - 2*beta**2 + 2j*beta*np.sqrt(1 - beta**2)

def layer_mats(h, lam, G, rho, mass='consistent'):
    """Return per-layer 2x2 blocks for x, z (in-plane) and y (antiplane) per Kausel Table 1.
    Ordering inside a layer: node1(top), node2(bottom)."""
    Mlin = np.array([[2, 1], [1, 2]]) / 6.0
    if mass == 'consistent':
        Mm = Mlin
    elif mass == 'lumped':
        Mm = np.eye(2) / 2.0
    elif mass == 'mixed':
        Mm = 0.5 * Mlin + 0.5 * np.eye(2) / 2.0
    Kd = np.array([[1, -1], [-1, 1]])
    lp2 = lam + 2*G
    Ax = h*lp2*Mlin; Az = h*G*Mlin; Ay = h*G*Mlin
    Gx = G/h*Kd; Gz = lp2/h*Kd; Gy = G/h*Kd
    M = rho*h*Mm
    # B_xz (rows x1,x2 ; cols z1,z2)
    Bxz = 0.5*np.array([[lam - G, -(lam + G)], [lam + G, -(lam - G)]])
    return Ax, Az, Ay, Gx, Gz, Gy, M, Bxz

def assemble(layers, mass='consistent', base='rigid', dashpot=None, omega=None):
    """layers: list of (h, Vs, Vp, rho, beta_s, beta_p). N interfaces = nl+1 (rigid base: last fixed).
    Returns global matrices with free dofs (for rigid base drop bottom interface)."""
    nl = len(layers); ni = nl + 1
    Ax = np.zeros((ni, ni), complex); Az = Ax.copy(); Ay = Ax.copy()
    Gx = Ax.copy(); Gz = Ax.copy(); Gy = Ax.copy(); M = Ax.copy(); Bxz = Ax.copy()
    for j, (h, Vs, Vp, rho, bs, bp) in enumerate(layers):
        Gs = rho*Vs**2*cmod(bs); Mp = rho*Vp**2*cmod(bp); lam = Mp - 2*Gs
        ax, az, ay, gx, gz, gy, m, bxz = layer_mats(h, lam, Gs, rho, mass)
        s = slice(j, j+2)
        Ax[s, s] += ax; Az[s, s] += az; Ay[s, s] += ay
        Gx[s, s] += gx; Gz[s, s] += gz; Gy[s, s] += gy; M[s, s] += m; Bxz[s, s] += bxz
    C_extra = {}
    if base == 'rigid':
        keep = slice(0, ni-1)
        mats = [X[keep, keep] for X in (Ax, Az, Ay, Gx, Gz, Gy, M, Bxz)]
        return mats
    elif base == 'dashpot':
        # Lysmer-Kuhlemeyer dashpots at bottom interface: add i*w*rho*V to diag of G (C=G-w^2M+iwD)
        rho_b, Vs_b, Vp_b = dashpot
        Gx[-1, -1] += 1j*omega*rho_b*Vs_b
        Gy[-1, -1] += 1j*omega*rho_b*Vs_b
        Gz[-1, -1] += 1j*omega*rho_b*Vp_b
        return [Ax, Az, Ay, Gx, Gz, Gy, M, Bxz]

def modes(mats, omega):
    Ax, Az, Ay, Gx, Gz, Gy, M, Bxz = mats
    n = Ax.shape[0]
    Cx = Gx - omega**2*M; Cz = Gz - omega**2*M; Cy = Gy - omega**2*M
    # Eq.(16): [k^2 Ax + Cx, Bxz; k^2 Bxz^T, k^2 Az + Cz] {phi_x; k phi_z} = 0
    # -> (k^2 Abar + Cbar) Z = 0 with Abar=[[Ax,0],[Bxz^T,Az]], Cbar=[[Cx,Bxz],[0,Cz]]
    Abar = np.block([[Ax, np.zeros((n, n))], [Bxz.T, Az]])
    Cbar = np.block([[Cx, Bxz], [np.zeros((n, n)), Cz]])
    w, Z = eig(-Cbar, Abar)   # w = k^2
    k = np.sqrt(w.astype(complex))
    # choose Im(k)<0, or Re(k)>0 if real
    k = np.where(np.imag(k) > 0, -k, k)
    k = np.where((np.abs(np.imag(k)) < 1e-14*np.abs(k)) & (np.real(k) < 0), -k, k)
    phx = Z[:n, :]; phz = Z[:, :][n:, :] / k
    # normalization Eq.(22a): Y^T Abar Z = k  with Y={k phx; phz}
    for j in range(2*n):
        Y = np.concatenate([k[j]*phx[:, j], phz[:, j]])
        Zj = np.concatenate([phx[:, j], k[j]*phz[:, j]])
        nrm = Y @ Abar @ Zj / k[j]
        s = 1/np.sqrt(nrm)
        phx[:, j] *= s; phz[:, j] *= s
    # Love
    wl, Vl = eig(-Cy, Ay)
    kl = np.sqrt(wl.astype(complex)); kl = np.where(np.imag(kl) > 0, -kl, kl)
    for j in range(n):
        nrm = Vl[:, j] @ Ay @ Vl[:, j]
        Vl[:, j] /= np.sqrt(nrm)
    return k, phx, phz, kl, Vl

def direct_flex(mats, omega, k):
    Ax, Az, Ay, Gx, Gz, Gy, M, Bxz = mats
    n = Ax.shape[0]
    K = np.block([[Ax*k**2 + Gx - omega**2*M, Bxz*k], [Bxz.T*k, Az*k**2 + Gz - omega**2*M]])
    return np.linalg.inv(K)

def modal_flex(md, k):
    kr, phx, phz, kl, phy = md
    D = 1/(k**2 - kr**2)
    Fxx = (phx*D) @ phx.T
    Fxz = k*(phx*(D/kr)) @ phz.T
    Fzx = k*(phz*(D/kr)) @ phx.T
    Fzz = (phz*D) @ phz.T
    return np.block([[Fxx, Fxz], [Fzx, Fzz]])

if __name__ == '__main__':
    import warnings; warnings.filterwarnings('ignore')
    np.set_printoptions(precision=5, linewidth=150)
    # ---- test 1: modal vs direct flexibility
    layers = [(1.0, 200, 400, 2.0, 0.05, 0.05)]*5 + [(2.0, 400, 800, 2.1, 0.03, 0.03)]*5
    om = 2*np.pi*10.0
    mats = assemble(layers, mass='consistent')
    md = modes(mats, om)
    for kk in [0.3, 1.7+0.2j]:
        Fd = direct_flex(mats, om, kk); Fm = modal_flex(md, kk)
        print('test1 modal vs direct rel err (k=%s):' % kk, np.abs(Fd - Fm).max()/np.abs(Fd).max())
    # Love modal vs direct
    Ax, Az, Ay, Gx, Gz, Gy, M, Bxz = mats
    kr, phx, phz, kl, phy = md
    kk = 0.9
    Fy_d = np.linalg.inv(Ay*kk**2 + Gy - om**2*M)
    Fy_m = (phy*(1/(kk**2 - kl**2))) @ phy.T
    print('test1b Love modal vs direct:', np.abs(Fy_d - Fy_m).max()/np.abs(Fy_d).max())
