"""Independent high-precision SG references (requires mpmath).

Uses the closed p=2 kernel in rapidity space, not production Fourier quadrature.
The pair is asymptotic Bethe-Yang, not an exact finite-volume excitation.
"""
import mpmath as mp

mp.mp.dps = 90


def phase(theta):
    return mp.quad(lambda x: x / (mp.pi * mp.sinh(x)) if x else 1 / mp.pi, [0, theta])


theta = mp.findroot(lambda x: 10 * mp.sinh(x) + phase(2 * x) - mp.pi, mp.mpf("0.3"))
print("chi(0.7,p=2) =", phase(mp.mpf("0.7")))
print("B1 energy(k=0.7,M=1,p=0.4) =", mp.hypot(2 * mp.sin(mp.pi * mp.mpf("0.4") / 2), mp.mpf("0.7")))
print("pair rapidity(M=1,L=10,p=2,I=1/2) =", theta)
print("pair energy =", 2 * mp.cosh(theta))
