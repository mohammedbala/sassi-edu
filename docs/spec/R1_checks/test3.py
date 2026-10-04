from greens import *
# graded profile, uniform material, deep rigid base
Vs=100.; nu=0.3; Vp=Vs*np.sqrt(2*(1-nu)/(1-2*nu)); rho=2.0; Gs=rho*Vs**2
hs=[]; z=0; h=0.05
while z<400:
    hs.append(h); z+=h; h=min(h*1.06,10.0)
layers=[(hh,Vs,Vp,rho,0.005,0.005) for hh in hs]
mats=assemble(layers,mass='consistent')
om=2*np.pi*0.01   # quasi-static (rigid base at 400 m, f1 = Vs/4H = 0.0625 Hz)
md=modes(mats,om)
print('N layers',len(hs))
for r in [1.0,2.0,5.0]:
    ur_h,ut_h,uz_h,ur_v,uz_v=green_point(md,0,0,r)
    # Boussinesq/Cerruti (halfspace surface, z down positive):
    # vertical load P down: uz_down=P(1-nu)/(2 pi G r); ur=-(1-2nu)P/(4 pi G r)
    # horizontal P_x: ux(theta=0)=P/(2piGr) ; ux(theta=90)=P(1-nu)/(2piGr); uz_down(theta=0) = -(1-2nu)P/(4 pi G r)
    B=1/(2*np.pi*Gs*r)
    print('r=%g'%r)
    print('  horiz load: u_rho~ %.4f (Cerruti %.4f), v~ %.4f (Cerruti %.4f), u_z_up~ %.4f (Cerruti up %.4f)'%(
        ur_h.real/B, 1.0, ut_h.real/B, (1-nu), uz_h.real/B, (1-2*nu)/2))
    # vertical load DOWN = p_z=-1 in z-up convention
    print('  vert load (down): u_rho %.4f (Bous %.4f), u_z_up %.4f (Bous up %.4f)'%(
        (-ur_v).real/B, -(1-2*nu)/2, (-uz_v).real/B, -(1-nu)))
