"""Reconstruction of SASSI POINT3: axisymmetric FE central zone (radius R0, one element radially)
+ cylindrical consistent transmitting boundary (Kausel/Waas) built from TLM eigen-solution.
Verifies far field against Kausel (1981) point-load Green functions (Eq. 88-89).
Kausel convention: z up, interfaces numbered from top, e^{iwt}, outgoing H^(2).
"""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from scipy.special import hankel2
from tlm import cmod, assemble, modes
from greens import green_point

def H(n, x): return hankel2(n, x)
def dH(n, x): return 0.5*(hankel2(n-1, x) - hankel2(n+1, x))
def ddH(n, x): return -dH(n, x)/x - (1 - n**2/x**2)*H(n, x)

def layer_props(layers):
    out = []
    for (h, Vs, Vp, rho, bs, bp) in layers:
        G = rho*Vs**2*cmod(bs); Mp = rho*Vp**2*cmod(bp); lam = Mp - 2*G
        out.append((h, lam, G, rho))
    return out

def EQ(layers, nfree, coef):
    """assembled E_c=int c N N dz and Q_c=int c N N' dz (z up, node1 top)."""
    ni = len(layers) + 1
    E = np.zeros((ni, ni), complex); Q = np.zeros((ni, ni), complex)
    for j, (h, lam, G, rho) in enumerate(layer_props(layers)):
        c = {'lp2': lam + 2*G, 'lam': lam, 'G': G}[coef]
        s = slice(j, j+2)
        E[s, s] += c*h/6*np.array([[2, 1], [1, 2]])
        Q[s, s] += c/2*np.array([[1, -1], [1, -1]])
    return E[:nfree, :nfree], Q[:nfree, :nfree]

def Psi(md, mu, r):
    kr, phx, phz, kl, phy = md
    xR = kr*r; xL = kl*r
    top = np.hstack([phx*dH(mu, xR), phy*(mu*H(mu, xL)/xL)])
    mid = np.hstack([phx*(mu*H(mu, xR)/xR), phy*dH(mu, xL)])
    bot = np.hstack([-phz*H(mu, xR), np.zeros_like(phy)])
    return np.vstack([top, mid, bot])

def dPsi(md, mu, r):
    kr, phx, phz, kl, phy = md
    xR = kr*r; xL = kl*r
    top = np.hstack([phx*(kr*ddH(mu, xR)), phy*(mu*kl*(dH(mu, xL)/xL - H(mu, xL)/xL**2))])
    mid = np.hstack([phx*(mu*kr*(dH(mu, xR)/xR - H(mu, xR)/xR**2)), phy*(kl*ddH(mu, xL))])
    bot = np.hstack([-phz*(kr*dH(mu, xR)), np.zeros_like(phy)])
    return np.vstack([top, mid, bot])

def transmitting_boundary(layers, md, mu, R0):
    n = md[1].shape[0]
    Elp, _ = EQ(layers, n, 'lp2'); El, Ql = EQ(layers, n, 'lam'); EG, QG = EQ(layers, n, 'G')
    P = Psi(md, mu, R0); dP = dPsi(md, mu, R0)
    ur, ut, uz = P[:n], P[n:2*n], P[2*n:]
    urr, utr, uzr = dP[:n], dP[n:2*n], dP[2*n:]
    fr = R0*(Elp@urr + El@(ur - mu*ut)/R0 + Ql@uz)
    ft = R0*(EG@(mu*ur/R0 + utr - ut/R0))
    fz = R0*(QG@ur + EG@uzr)
    T = np.vstack([fr, ft, fz])
    return -T@np.linalg.inv(P)

def core_matrices(layers, nfree, mu, R0, omega, mass='consistent', ngr=3):
    """Axisymmetric-harmonic FE core, one element radially (nodes rho=0 and rho=R0), linear in z.
    dof order: [ur(axis,all ifc), ut(axis), uz(axis), ur(R0), ut(R0), uz(R0)] -> later reduce.
    returns per-radian dynamic stiffness K - w^2 M (6n x 6n) in order [axis(r,t,z), boundary(r,t,z)]."""
    ni = len(layers) + 1
    nd = 6*ni
    K = np.zeros((nd, nd), complex); M = np.zeros((nd, nd), complex)
    g, wg = np.polynomial.legendre.leggauss(ngr)
    def idx(node_r, comp, ifc):  # node_r 0=axis,1=boundary; comp 0=r,1=t,2=z
        return (node_r*3 + comp)*ni + ifc
    for j, (h, lam, G, rho) in enumerate(layer_props(layers)):
        Dm = np.zeros((6, 6), complex)
        Dm[:3, :3] = lam; Dm[0, 0] = Dm[1, 1] = Dm[2, 2] = lam + 2*G
        Dm[3, 3] = Dm[4, 4] = Dm[5, 5] = G
        for gr, wr in zip(g, wg):
            r = 0.5*R0*(gr + 1); wr_ = 0.5*R0*wr
            L = np.array([1 - r/R0, r/R0]); dL = np.array([-1/R0, 1/R0])
            for gz, wz in zip(g, wg):
                s = 0.5*(gz + 1)       # 0 at top node, 1 at bottom node
                Nz = np.array([1 - s, s]); dNz = np.array([1/h, -1/h])  # d/dz, z up: top node value increases upward
                wz_ = 0.5*h*wz
                B = np.zeros((6, nd), complex); Nmat = np.zeros((3, nd))
                for a in range(2):        # radial node
                    for b in range(2):    # z node
                        ifc = j + b
                        Nr = L[a]*Nz[b]; Nr_r = dL[a]*Nz[b]; Nr_z = L[a]*dNz[b]
                        ir, it, iz = idx(a, 0, ifc), idx(a, 1, ifc), idx(a, 2, ifc)
                        # eps = [e_rr, e_tt, e_zz, g_rz, g_rt~, g_tz~]
                        B[0, ir] += Nr_r
                        B[1, ir] += Nr/r; B[1, it] += -mu*Nr/r
                        B[2, iz] += Nr_z
                        B[3, ir] += Nr_z; B[3, iz] += Nr_r
                        B[4, ir] += mu*Nr/r; B[4, it] += Nr_r - Nr/r
                        B[5, it] += Nr_z; B[5, iz] += mu*Nr/r
                        Nmat[0, ir] += Nr; Nmat[1, it] += Nr; Nmat[2, iz] += Nr
                wgt = wr_*wz_*r
                K += B.T@Dm@B*wgt
                M += rho*(Nmat.T@Nmat)*wgt
    if mass == 'lumped_z':
        pass
    return K - omega**2*M, ni

def solve_point(layers, md, mu, R0, omega, load_ifc, load_dir):
    """Unit point load at axis node of interface load_ifc. load_dir: 'x' (mu=1) or 'z' (mu=0).
    rigid base: bottom interface fixed. returns axis displacement vector and boundary dofs."""
    n = md[1].shape[0]          # free interfaces
    Kc, ni = core_matrices(layers, n, mu, R0, omega)
    # dof selection: free interfaces only (exclude bottom ifc if rigid)
    def idx(node_r, comp, ifc): return (node_r*3 + comp)*ni + ifc
    free_ifc = range(n)
    # build reduced dof list with axis constraints
    # unknowns: axis: mu=0 -> uz ; mu=1 -> ur(=ut)  ; boundary: ur,ut,uz
    cols = []   # list of lists of (full index, coefficient)
    for i in free_ifc:
        if mu == 0:
            cols.append([(idx(0, 2, i), 1.0)])
        else:
            cols.append([(idx(0, 0, i), 1.0), (idx(0, 1, i), 1.0)])
    nax = len(cols)
    for comp in range(3):
        for i in free_ifc:
            cols.append([(idx(1, comp, i), 1.0)])
    T = np.zeros((Kc.shape[0], len(cols)))
    for c, lst in enumerate(cols):
        for (ii, v) in lst: T[ii, c] = v
    Kr = T.T@Kc@T
    Rtb = transmitting_boundary(layers, md, mu, R0)
    Kr[nax:, nax:] += Rtb
    F = np.zeros(len(cols), complex)
    F[load_ifc] = 1.0/(np.pi if mu >= 1 else 2*np.pi)
    U = np.linalg.solve(Kr, F)
    Uax = U[:nax]; Ub = U[nax:]
    alpha = np.linalg.solve(Psi(md, mu, R0), Ub)
    return Uax, Ub, alpha

if __name__ == '__main__':
    layers = [(0.5, 150, 300, 1.9, 0.05, 0.05)]*8 + [(1.0, 250, 500, 2.0, 0.04, 0.04)]*10 + [(2.0, 400, 800, 2.1, 0.03, 0.03)]*10
    f = 5.0; om = 2*np.pi*f
    mats = assemble(layers, mass='consistent')
    md = modes(mats, om)
    n = md[1].shape[0]
    R0 = 0.9
    for (mu, d) in [(1, 'x'), (0, 'z')]:
        for load_ifc in [0, 6]:
            Uax, Ub, alpha = solve_point(layers, md, mu, R0, om, load_ifc, d)
            print('mu=%d load at ifc %d' % (mu, load_ifc))
            for r in [1.8, 3.0, 6.0, 12.0, 25.0]:
                P = Psi(md, mu, r)@alpha
                ur, ut, uz = P[:n], P[n:2*n], P[2*n:]
                for m in [0, 4]:
                    g = green_point(md, m, load_ifc, r)
                    if mu == 1:
                        ref = (g[0], g[1], g[2]); got = (ur[m], ut[m], uz[m])
                    else:
                        ref = (g[3], 0, g[4]); got = (ur[m], ut[m], uz[m])
                    rel = [abs(a - b)/max(abs(b), 1e-30) for a, b in zip(got, ref) if abs(b) > 0]
                    print('   r=%5.1f obs ifc %d: FEcore+TB vs Kausel point: rel err %s' % (r, m, np.round(rel, 4)))
            print('   axis displacement at load ifc (self-flexibility):', Uax[load_ifc])
