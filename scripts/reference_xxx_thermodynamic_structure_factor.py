#!/usr/bin/env python3
"""Optional mpmath oracle: different tail subtraction and tanh-sinh meshes.

Print I(rho) for the C++ regression points. Not a runtime dependency.
The C++ code subtracts 2*(1-exp(-2*x))*cos(a*x)/x; here use exp(-x),
whose analytic integral is log(1+1/a^2), and integrate out to 160.
"""
import mpmath as mp

def transition_integral(rho):
    a = 4 * rho

    def remainder(x):
        t = mp.tanh(x)
        return ((1+t)/mp.tanh(2*x)*(-2*mp.sin(a*x/2)**2)
                + t+t*t + 2*mp.expm1(-x)*mp.cos(a*x))/x

    return mp.log1p(1/a**2) + mp.quad(remainder, [mp.mpf(i)/2 for i in range(321)])


mp.mp.dps = 65
for text in ('0.0001', '0.1', '0.5', '1', '2', '5'):
    print(text, mp.nstr(transition_integral(mp.mpf(text)), 52), flush=True)
q, omega = mp.mpf('1.5'), mp.mpf('1.8')
lower, upper = mp.pi/2*mp.sin(q), mp.pi*mp.sin(q/2)
rho = mp.acosh(mp.sqrt((upper**2-lower**2)/(omega**2-lower**2)))/mp.pi
density = mp.exp(-transition_integral(rho))/(2*mp.sqrt(upper**2-omega**2))
print('Szz(1.5,1.8)', mp.nstr(density, 55))
