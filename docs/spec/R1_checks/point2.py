"""Reconstruction of SASSI POINT2 (plane strain, in-plane P-SV): FE central strip -R0..R0 (2 elements wide,
linear in z per TLM layer) + Waas transmitting boundaries at x=+-R0. Verified vs Kausel (1981) Eq.(47)."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import assemble, modes, cmod
from point3 import layer_props, EQ

def tb_plane(layers, md, side):
    """Return physical-variable TB stiffness (2n x 2n), dof order [ux(ifc..), uz(ifc..)], z up.
    side=+1 right boundary (outgoing e^{-ik(x-xb)}), -1 left boundary (outgoing e^{+ik(x-xb)})."""
    kr, phx, phz, kl, phy = md
    n = phx.shape[0]
    Elp, _ = EQ(layers, n, 'lp2'); El, Ql = EQ(layers, n, 'lam'); EG, QG = EQ(layers, n, 'G')
    # physical mode shapes at boundary: ux = phx, uz = -i*s*phz (s=+1 right-going; -1 left-going)
    s = side
    ux = phx.astype(complex); uz = -1j*s*phz
    # x-derivatives: d/dx -> -i*s*k
    dux = -1j*s*kr*ux; duz = -1j*s*kr*uz
    sxx = Elp@dux + Ql@uz          # int N^T sigma_xx dz
    txz = QG@ux + EG@duz           # int N^T tau_xz dz
    T = s*np.vstack([sxx, txz])    # force on interior (outward normal = s*e_x)
    Psi = np.vstack([ux, uz])
    return -T@np.linalg.inv(Psi), Psi

def core_plane(layers, nfree, R0, omega, ngr=3):
    ni = len(layers) + 1
    xs = [-R0, 0.0, R0]
    nn = 3*ni
    nd = 2*nn
    K = np.zeros((nd, nd), complex); M = np.zeros((nd, nd), complex)
    def idx(ix, comp, ifc): return comp*nn + ix*ni + ifc
    g, wg = np.polynomial.legendre.leggauss(ngr)
    for j, (h, lam, G, rho) in enumerate(layer_props(layers)):
        D = np.array([[lam + 2*G, lam, 0], [lam, lam + 2*G, 0], [0, 0, G]])
        for e in range(2):
            x0, x1 = xs[e], xs[e+1]; le = x1 - x0
            for gx, wx in zip(g, wg):
                xi = 0.5*(gx + 1); L = np.array([1 - xi, xi]); dL = np.array([-1/le, 1/le])
                for gz, wz in zip(g, wg):
                    s = 0.5*(gz + 1); Nz = np.array([1 - s, s]); dNz = np.array([1/h, -1/h])
                    B = np.zeros((3, nd), complex); Nm = np.zeros((2, nd))
                    for a in range(2):
                        for b in range(2):
                            ix = e + a; ifc = j + b
                            iu, iw = idx(ix, 0, ifc), idx(ix, 1, ifc)
                            B[0, iu] += dL[a]*Nz[b]; B[1, iw] += L[a]*dNz[b]
                            B[2, iu] += L[a]*dNz[b]; B[2, iw] += dL[a]*Nz[b]
                            Nm[0, iu] += L[a]*Nz[b]; Nm[1, iw] += L[a]*Nz[b]
                    wgt = 0.5*le*wx*0.5*h*wz
                    K += B.T@D@B*wgt; M += rho*Nm.T@Nm*wgt
    return K - omega**2*M, ni, idx

if __name__ == '__main__':
    layers = [(1.0, 150, 300, 1.9, 0.05, 0.05)]*6 + [(2.0, 300, 600, 2.0, 0.05, 0.05)]*6
    f = 4.0; om = 2*np.pi*f
    md = modes(assemble(layers), om); kr, phx, phz, kl, phy = md; n = phx.shape[0]
    R0 = 1.0
    Kc, ni, idx = core_plane(layers, n, R0, om)
    keep = [idx(ix, c, i) for c in range(2) for ix in range(3) for i in range(n)]
    Kr = Kc[np.ix_(keep, keep)]
    pos = {k: p for p, k in enumerate(keep)}
    Rr, Pr = tb_plane(layers, md, +1); Rl, Pl = tb_plane(layers, md, -1)
    def bnd(ix): return [pos[idx(ix, c, i)] for c in range(2) for i in range(n)]
    Kr[np.ix_(bnd(2), bnd(2))] += Rr; Kr[np.ix_(bnd(0), bnd(0))] += Rl
    for comp, name in [(0, 'x'), (1, 'z')]:
        F = np.zeros(len(keep), complex); load_ifc = 3
        F[pos[idx(1, comp, load_ifc)]] = 1.0
        U = np.linalg.solve(Kr, F)
        Ub = U[bnd(2)]; a = np.linalg.solve(Pr, Ub)
        for x in [2.0, 5.0, 12.0]:
            Ux = (phx*np.exp(-1j*kr*(x - R0)))@a; Uz = (-1j*phz*np.exp(-1j*kr*(x - R0)))@a
            m = 0
            # Kausel Eq.(47) physical: u_x=(1/2i)sum phx phx E/k (x-load); vertical comps carry i factors
            E = np.exp(-1j*kr*abs(x))
            if comp == 0:
                ref_x = np.sum(phx[m]*phx[load_ifc]*E/kr)/(2j)
                ref_z = -1j*np.sum(phz[m]*phx[load_ifc]*E/kr)/(2j)   # i*uz = scaled -> uz = -i*scaled
            else:
                ref_x = 1j*np.sum(phx[m]*phz[load_ifc]*E/kr)/(2j)     # scaled load i*pz -> factor i
                ref_z = np.sum(phz[m]*phz[load_ifc]*E/kr)/(2j)
            print('%s-load x=%5.1f: ux FE+TB %s vs Eq47 %s | uz %s vs %s' % (name, x, np.round(Ux[m], 9), np.round(ref_x, 9), np.round(Uz[m], 9), np.round(ref_z, 9)))
