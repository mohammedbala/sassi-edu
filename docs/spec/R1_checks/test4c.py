import numpy as np, warnings
warnings.filterwarnings('ignore')
from tlm import cmod
from test4 import prof, hs, shake_tf
from test4b import column
from scipy.optimize import brentq
def tlm_outcrop2(f,nhs,frac=0.05,mass='consistent',dist='uniform',total=1.5):
    w=2*np.pi*f; hl=[];props=[]
    for h,Vs,rho,b in prof:
        n=max(1,int(np.ceil(h/(frac*Vs/f)))); hl+=[h/n]*n; props+=[(Vs,rho,b,mass)]*n
    Vh,rh,bh=hs; D=total*Vh/f
    if dist=='uniform': hh=[D/nhs]*nhs
    else:
        h1=min(hl[-1]*Vh/prof[-1][1], D/nhs)  # first HS sublayer scaled by velocity ratio
        q=brentq(lambda q: h1*(q**nhs-1)/(q-1)-D,1.0000001,10); hh=[h1*q**i for i in range(nhs)]
    c=rh*Vh*np.sqrt(cmod(bh))
    U1=column(hl+hh,props+[(Vh,rh,bh,mass)]*nhs,w,c); U2=column(hh,[(Vh,rh,bh,mass)]*nhs,w,c)
    return U1[0]/U2[0]
for dist in ['uniform','geom_scaled']:
  print(dist)
  for f in [0.5,1.7,4.,10.,15.]:
    s=shake_tf(f); row=[]
    for nhs in [5,10,20]:
        t=tlm_outcrop2(f,nhs,dist=dist); row.append('n=%2d %.4f/%5.1fdeg'%(nhs,abs(t)/abs(s),np.degrees(np.angle(t/s))))
    print('  f=%5.2f '%f+'  '.join(row))
