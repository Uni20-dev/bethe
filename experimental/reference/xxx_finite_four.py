#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Independent full-root N16/N64 finite-four form-factor controls.

Requires mpmath >= 1.2. Uses direct multiplicative equations and unreduced
complex determinants; no C++ stabilized matrices or private data are imported.
The N64 deviation requires more than 100 digits merely to resolve its roots.
"""
import argparse
from itertools import combinations
import json
import mpmath as mp
from xxx_ground_fixture import ground


def bethe_residual(n, x, others):
    return ((x+mp.mpc(0, '.5'))/(x-mp.mpc(0, '.5')))**n \
        * mp.fprod((x-y-1j)/(x-y+1j) for y in others)-1


def max_residual(n, roots):
    return max(abs(bethe_residual(n, x, roots[:j]+roots[j+1:]))
               for j, x in enumerate(roots))


def energy(n, roots):
    value = n/mp.mpf(4)-sum(mp.mpf('.5')/(x*x+mp.mpf('.25')) for x in roots)
    if abs(mp.im(value)) > mp.mpf('1e-30'):
        raise ValueError('nonreal energy')
    return mp.re(value)


def takahashi_residual(n, labels, branch, roots, v):
    def theta(k,x):
        return 2*mp.atan(2*x/k)
    out=[n*theta(1,a)-2*mp.pi*i-sum(theta(2,a-b) for b in roots)
         -theta(1,a-v)-theta(3,a-v) for a,i in zip(roots,labels)]
    return out+[n*theta(2,v)-2*mp.pi*branch
                -sum(theta(1,v-a)+theta(3,v-a) for a in roots)]


def ideal_string_roots(n, holes, branch):
    if n < 8 or n % 4 or branch not in (-1, 0, 1):
        raise ValueError('require N>=8 divisible by four and branch -1,0,1')
    if (len(holes) != 4 or len(set(holes)) != 4 or
            any(int(x) != x or abs(x) > n//4 for x in holes)):
        raise ValueError('require four distinct integer holes in [-N/4,N/4]')
    if mp.mp.dps < 50:
        raise ValueError('use at least 50 decimal digits, and more for long strings')
    labels = [i for i in range(-n//4, n//4+1) if i not in holes]

    def takahashi(*x):
        return takahashi_residual(n,labels,branch,x[:-1],x[-1])

    guess = [mp.tan(mp.pi*i/n)/2 for i in labels]+[mp.mpf(branch)/2]
    seed = list(mp.findroot(takahashi, guess,
                           J=lambda *x: string_gaudin_matrix(n,list(x[:-1]),x[-1]),
                           tol=mp.mpf('1e-45'), maxsteps=100))
    if max(abs(x) for x in takahashi(*seed)) > mp.mpf('1e-30'):
        raise ValueError('Takahashi roots unconverged')
    return seed


def solve(n, holes, branch):
    seed = ideal_string_roots(n,holes,branch)
    reals, v = seed[:-1], seed[-1]
    if abs(v) < mp.mpf('1e-25'):
        raise ValueError('singular centred pair is outside this research solver')
    ideal = reals+[v+mp.mpc(0,'.5'), v-mp.mpc(0,'.5')]
    r = bethe_residual(n, ideal[-2], reals)+1
    deviation = mp.re(1/(r-1))

    solution = list(mp.findroot(lambda *x: deviated_equations(n,x), reals+[v,deviation],
                               J=lambda *x: deviated_jacobian(n,x),
                               tol=mp.mpf('1e-42'), maxsteps=100))
    centre, delta = solution[-2:]
    if mp.mpf('.5')+delta <= 0:
        raise ValueError('collapsed or reversed string is not supported')
    upper = centre+1j*(mp.mpf('.5')+delta)
    roots = solution[:-2]+[upper,mp.conj(upper)]
    error = max_residual(n, roots)
    if error > mp.mpf('1e-30'):
        raise ValueError('full complex Bethe equations failed; increase precision or reject state')
    phase = mp.fprod((x+mp.mpc(0,'.5'))/(x-mp.mpc(0,'.5')) for x in roots)
    momentum_index = (sum(holes)-branch) % n
    if abs(phase-mp.exp(2j*mp.pi*momentum_index/n)) > mp.mpf('1e-28'):
        raise ValueError('solution left its expected momentum sector')
    return dict(roots=roots, ideal_roots=ideal, centre=centre, deviation=delta,
                residual=error, energy=energy(n,roots), momentum_index=momentum_index)


def deviated_equations(n, variables):
    roots,centre,delta=variables[:-2],variables[-2],variables[-1]
    upper=centre+1j*(mp.mpf('.5')+delta)
    lower=mp.conj(upper)
    out=[mp.im(bethe_residual(n,a,list(roots[:j])+list(roots[j+1:])+[upper,lower]))
         for j,a in enumerate(roots)]
    pair=(bethe_residual(n,upper,roots)+1)*delta/(1+delta)-1
    return out+[mp.re(pair),mp.im(pair)]


def deviated_jacobian(n, variables):
    """Exact real Jacobian, preserving delta/(1+delta) without subtraction."""
    roots,centre,delta=variables[:-2],variables[-2],variables[-1]
    upper=centre+1j*(mp.mpf('.5')+delta)
    lower=mp.conj(upper)
    size=len(roots)
    jac=mp.matrix(size+2)
    for i,x in enumerate(roots):
        phase=bethe_residual(n,x,list(roots[:i])+list(roots[i+1:])+[upper,lower])+1
        kp,km=2/(1+(x-upper)**2),2/(1+(x-lower)**2)
        diagonal=n/(x*x+mp.mpf('.25'))-kp-km
        for j,y in enumerate(roots):
            if i!=j:
                kernel=2/(1+(x-y)**2)
                diagonal-=kernel
                jac[i,j]=mp.im(-1j*phase*kernel)
        jac[i,i]=mp.im(-1j*phase*diagonal)
        jac[i,size]=mp.im(-1j*phase*(kp+km))
        jac[i,size+1]=mp.im(phase*(kp-km))
    phase=(bethe_residual(n,upper,roots)+1)*delta/(1+delta)
    kernels=[2/(1+(upper-x)**2) for x in roots]
    diagonal=n/(upper*upper+mp.mpf('.25'))-sum(kernels)
    derivatives=[-1j*phase*k for k in kernels]
    derivatives += [-1j*phase*diagonal,phase*(diagonal+1/(delta*(1+delta)))]
    for j,z in enumerate(derivatives):
        jac[size,j],jac[size+1,j]=mp.re(z),mp.im(z)
    return jac


def absolute_determinant(matrix):
    """Equilibrate columns before mpmath's norm-dependent pivot threshold."""
    scaled = matrix.copy()
    factors = []
    for j in range(matrix.cols):
        scale = max(abs(matrix[i,j]) for i in range(matrix.rows))
        if scale == 0:
            raise ValueError('zero determinant column; weight unresolved')
        factors.append(scale)
        for i in range(matrix.rows):
            scaled[i,j] /= scale
    result = abs(mp.det(scaled))*mp.fprod(factors)
    if not mp.isfinite(result) or result == 0:
        raise ValueError('unresolved determinant; increase precision')
    return result


def gaudin_matrix(n, roots):
    matrix = mp.matrix(len(roots))
    for i,x in enumerate(roots):
        for j,y in enumerate(roots):
            if i != j:
                matrix[i,j] = 2/(1+(x-y)**2)
        matrix[i,i] = n/(x*x+mp.mpf('.25'))-sum(matrix[i,j] for j in range(len(roots)) if j!=i)
    return matrix


def gaudin_determinant(n, roots):
    return absolute_determinant(gaudin_matrix(n,roots))


def string_gaudin_matrix(n, reals, centre):
    """Jacobian of the ideal one-/two-string Takahashi equations."""
    def kernel(k,x):
        return k/(x*x+mp.mpf(k*k)/4)
    size=len(reals)
    matrix=mp.matrix(size+1)
    for i,x in enumerate(reals):
        for j,y in enumerate(reals):
            if i!=j:
                matrix[i,j]=kernel(2,x-y)
        matrix[i,size]=matrix[size,i]=kernel(1,x-centre)+kernel(3,x-centre)
        matrix[i,i]=n*kernel(1,x)-sum(matrix[i,j] for j in range(size+1) if j!=i)
    matrix[size,size]=n*kernel(2,centre)-sum(matrix[size,j] for j in range(size))
    return matrix


def raising_matrix(n, mu, lam):
    """Fill the supplied rapidity columns and the local-operator column."""
    h=mp.matrix(len(mu))
    for b,v in enumerate(lam):
        phase=((v+mp.mpc(0,'.5'))/(v-mp.mpc(0,'.5')))**n
        minus=mp.fprod(x-v-1j for x in mu)
        plus=mp.fprod(x-v+1j for x in mu)
        for a,u in enumerate(mu):
            h[a,b]=(minus/(u-v-1j)-phase*plus/(u-v+1j))/(u-v)
    for a,u in enumerate(mu):
        h[a,len(mu)-1]=1/(u*u+mp.mpf('.25'))
    return h


def raising_weight(n, mu, lam):
    """N times normalized local raising weight, unbroadened, without 2*pi.

    Full unreduced determinant with resolved conjugate roots. Unlike the real
    specialization, the plus product is NOT the conjugate of the minus product.
    """
    if len(mu) != n//2 or len(lam) != n//2-1:
        raise ValueError('wrong root counts')
    if max(max_residual(n,mu),max_residual(n,lam)) > mp.mpf('1e-30'):
        raise ValueError('form factor requires converged on-shell roots')
    h = raising_matrix(n,mu,lam)
    value = n*absolute_determinant(h)**2/(gaudin_determinant(n,mu)*gaudin_determinant(n,lam))
    value *= mp.fprod(abs(x*x+mp.mpf('.25')) for x in mu)
    value /= mp.fprod(abs(x*x+mp.mpf('.25')) for x in lam)
    for roots in (mu,lam):
        value /= mp.fprod(abs(1+(x-y)**2) for x,y in combinations(roots,2))
    if not mp.isfinite(value) or not 0 < value <= n:
        raise ValueError('invalid form-factor weight')
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits',type=int,default=150)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.digits<150:
        parser.error('use at least 150 decimal digits for the narrow N64 pair')
    cases=[(16,[-1,0,2,3],0,'-7.56318282038341941997737176833694899404308444e-12',
            '8.28523087152927669547401134062756665438561737e-8'),
           (64,[-14,-12,-11,-10],1,'1.41905199057999036182880353262087001032346067e-94',
            '1.33985918781507196786052358985966929800159504e-11')]
    with mp.workdps(args.digits):
        for n,holes,branch,expected_delta,expected_weight in cases:
            mu=ground(n)
            state=solve(n,holes,branch)
            weight=raising_weight(n,mu,state['roots'])
            if args.verify:
                for value,expected in [(state['deviation'],expected_delta),(weight,expected_weight)]:
                    if abs(value/mp.mpf(expected)-1)>mp.mpf('1e-40'):
                        raise RuntimeError(f'N{n} fixture differs from original-equation reference')
            print(json.dumps({'sites':n,'holes':holes,'branch':branch,'digits':args.digits,
                              'deviation':mp.nstr(state['deviation'],60),
                              'raising_weight':mp.nstr(weight,60),
                              'energy':mp.nstr(state['energy'],60),
                              'rational_residual':mp.nstr(state['residual'],10)}),flush=True)


if __name__=='__main__':
    main()
