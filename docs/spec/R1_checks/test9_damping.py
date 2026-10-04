"""Damping-convention check against Ostadan, Deng & Roesset (2004) Table 1 (fixed-base SDOF run in SASSI2000):
given beta 5/10/15/20 %  ->  '1/(2*Umax)' column 5.0/9.9/14.7/19.5 %, 'Half band' 4.8/9.7/15.9/22.3 %."""
import numpy as np
w = np.linspace(0.5, 1.5, 400001)
print('beta  form                    1/(2Hmax)   half-band')
for b in [0.05, 0.10, 0.15, 0.20]:
    for name, ks in [('G(1+2ib)', 1 + 2j*b), ('G(1-2b^2+2ib*sqrt(1-b^2))', 1 - 2*b*b + 2j*b*np.sqrt(1 - b*b))]:
        H = ks/(ks - w**2); Hm = abs(H).max()
        idx = np.where(abs(H) >= Hm/np.sqrt(2))[0]
        bw = (w[idx[-1]] - w[idx[0]])/(2*w[abs(H).argmax()])
        print('%.2f  %-26s %.4f     %.4f' % (b, name, 1/(2*Hm), bw))
