# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Selected nonsingular six-spinon triplets; arbitrary-precision test oracle.

Two two-strings or one three-string, with resolved finite-size deviations.
Requires mpmath >= 1.2; no private data or Bethe extension is imported.
Use --verify to check all eight committed root/energy/weight fixtures against fresh solves.
Includes the two N=12 non-self-conjugate-string configurations identified in
Deguchi-Giri (2015), Table I; not a general exceptional-string solver.
Not a thermodynamic spectral function. All energies
use H=sum S.S; weights are N|<e|S_0^+|0>|^2, without 2*pi.
"""
import argparse
import json

import mpmath as mp

from xxx_regular_roots import equations as original_equations, gaudin as gaudin_matrix, determinant


def bethe_residual(n, x, others):
    return ((x+.5j)/(x-.5j))**n * mp.fprod((x-y-1j)/(x-y+1j) for y in others)-1


def max_residual(n, roots):
    return max(abs(x) for x in original_equations(n, roots))


def energy(n, roots):
    value = n/mp.mpf(4)-sum(mp.mpf('.5')/(x*x+mp.mpf('.25')) for x in roots)
    if abs(mp.im(value)) > mp.mpf(10)**(-mp.mp.dps+25):
        raise ValueError('nonreal energy')
    return mp.re(value)


def ground(n):
    labels = [mp.mpf(j)-(mp.mpf(n)/2-1)/2 for j in range(n//2)]
    return list(mp.findroot(
        lambda *x: [2*n*mp.atan(2*a)-2*mp.pi*q-sum(2*mp.atan(a-b) for b in x)
                    for a,q in zip(x,labels)],
        [mp.tan(mp.pi*q/n)/2 for q in labels],
        J=lambda *x: gaudin_matrix(n,x), tol=mp.mpf(10)**(-mp.mp.dps+15), maxsteps=100))


def raising_weight(n, mu, roots):
    h=mp.matrix(len(mu))
    for b,v in enumerate(roots):
        phase=((v+.5j)/(v-.5j))**n
        minus=mp.fprod(u-v-1j for u in mu); plus=mp.fprod(u-v+1j for u in mu)
        for a,u in enumerate(mu):
            h[a,b]=(minus/(u-v-1j)-phase*plus/(u-v+1j))/(u-v)
    for a,u in enumerate(mu):
        h[a,len(mu)-1]=1/(u*u+mp.mpf('.25'))
    weight=n*determinant(h)**2/(determinant(gaudin_matrix(n,mu))*determinant(gaudin_matrix(n,roots)))
    weight*=mp.fprod(abs(u*u+mp.mpf('.25')) for u in mu)
    weight/=mp.fprod(abs(v*v+mp.mpf('.25')) for v in roots)
    for family in (mu,roots):
        weight/=mp.fprod(abs(1+(x-y)**2) for i,x in enumerate(family) for y in family[:i])
    return weight


def configuration(n, family, holes, branches):
    if n < 12 or n % 4 or family not in ('22', '3'):
        raise ValueError('require N>=12 divisible by four and family 22 or 3')
    holes, branches = sorted(map(mp.mpf, holes)), sorted(map(mp.mpf, branches))
    edge = mp.mpf(n)/4 + (mp.mpf('.5') if family == '3' else 0)
    window = [-edge+j for j in range(int(2*edge)+1)]
    allowed = list(map(mp.mpf, ('-1.5', '-.5', '.5', '1.5') if family == '22'
                       else ('-1', '0', '1')))
    if len(holes) != 6 or len(set(holes)) != 6 or any(h not in window for h in holes):
        raise ValueError('require six distinct holes in the family counting window')
    if (len(branches) != (2 if family == '22' else 1) or
            len(set(branches)) != len(branches) or any(b not in allowed for b in branches)):
        raise ValueError('invalid string quantum numbers')
    labels = [x for x in window if x not in holes] + branches
    lengths = [1]*(len(labels)-len(branches)) + ([2, 2] if family == '22' else [3])
    momentum = int(n/2-sum(labels)) % n
    return lengths, labels, momentum


def scattering(a, b, x, derivative=False):
    def theta(k):
        return k/(x*x+mp.mpf(k*k)/4) if derivative else 2*mp.atan(2*x/k)
    return ((theta(abs(a-b)) if a != b else 0) + theta(a+b)
            + 2*sum(theta(k) for k in range(abs(a-b)+2, a+b, 2)))


def takahashi(n, lengths, labels, centres, derivative=False):
    if derivative:
        g = mp.matrix(len(centres))
        for i, (a, x) in enumerate(zip(lengths, centres)):
            for j, (b, y) in enumerate(zip(lengths, centres)):
                if i != j:
                    g[i,j] = scattering(a, b, x-y, True)
            g[i,i] = n*a/(x*x+mp.mpf(a*a)/4)-sum(g[i,j] for j in range(len(centres)) if j != i)
        return g
    return [2*n*mp.atan(2*x/a)-2*mp.pi*label
            -sum(scattering(a,b,x-y) for j,(b,y) in enumerate(zip(lengths,centres)) if i != j)
            for i,(a,x,label) in enumerate(zip(lengths,centres,labels))]


def unpack(variables, real_count):
    roots = list(variables[:real_count])
    for i in range(real_count, len(variables), 2):
        u = variables[i]+1j*variables[i+1]
        roots += [u, mp.conj(u)]
    return roots


def resolved_equations(n, variables, real_count, derivative=False):
    roots = unpack(variables, real_count)
    phases = [bethe_residual(n,z,roots[:i]+roots[i+1:])+1 for i,z in enumerate(roots)]
    if not derivative:
        out = [mp.im(phases[i]) for i in range(real_count)]
        for i in range(real_count, len(roots), 2):
            out += [mp.re(phases[i])-1, mp.im(phases[i])]
        return out
    g = gaudin_matrix(n,roots)
    jac = mp.matrix(len(variables))
    for i in list(range(real_count))+list(range(real_count,len(roots),2)):
        derivatives = [-1j*phases[i]*g[i,j] for j in range(real_count)]
        for j in range(real_count,len(roots),2):
            derivatives += [-1j*phases[i]*(g[i,j]+g[i,j+1]), phases[i]*(g[i,j]-g[i,j+1])]
        for j,z in enumerate(derivatives):
            if i < real_count:
                jac[i,j] = mp.im(z)
            else:
                jac[i,j], jac[i+1,j] = mp.re(z), mp.im(z)
    return jac


def nonself_conjugate_sign(n, family, labels):
    """Only the two established Table-I configurations, not a generic fallback."""
    if n == 12 and family == '3':
        if labels == [mp.mpf('1.5'), mp.mpf('2.5'), mp.mpf(1)]:
            return 1
        if labels == [mp.mpf('-2.5'), mp.mpf('-1.5'), mp.mpf(-1)]:
            return -1
    return 0


def bethe_quantum_numbers(n, roots):
    return [mp.re((2*n*mp.atan(2*x)-sum(2*mp.atan(x-y) for j,y in enumerate(roots) if j != i))/(2*mp.pi))
            for i,x in enumerate(roots)]


def recovered_quantum_numbers(n, family, holes, branches, roots, *, unwrapped=False):
    """Recover labels, not just momentum. Non-self-conjugate mapping is local.

    Hagemans-Caux (2007), Eqs. 25/31 for individual self-conjugate strings;
    Deguchi-Giri (2015), Table I and Sec. VI for the two N=12 exceptions.
    The latter cannot be relabelled as the two-two-string family.
    """
    lengths,labels,_ = configuration(n,family,holes,branches)
    sign = nonself_conjugate_sign(n, family, labels)
    if sign:
        if len(roots) != 5:
            raise ValueError('non-self-conjugate configuration requires five roots')
        order = [0,1,2,3,4] if sign > 0 else [0,2,1,4,3]
        canonical = [sign*roots[i] for i in order]
        js = bethe_quantum_numbers(n, canonical)
        # Canonical order: real, inner upper/lower, outer upper/lower.
        # Table I assigns inner lower to the nominal one-string, and inner
        # upper to the three-string. The limiting label correction is -1/2
        # for each one-string and +1 for the three-string, modulo N.
        reals = sorted([sign*(js[0]-mp.mpf('.5')), sign*(js[2]-mp.mpf('.5'))])
        triple = sign*(js[1]+js[3]+js[4]-n+1)
        return reals+[triple]
    reals = lengths.count(1)
    js = bethe_quantum_numbers(n, roots)
    numbers = [js[i]+(mp.sign(roots[i]-roots[reals])/2 if family=='3' else 0) for i in range(reals)]
    if family == '22':
        numbers += [(1 if mp.im(roots[j])>mp.mpf('.5') else 0)+2*js[j]-mp.mpf(n)/2
                    +mp.sign(mp.re(roots[j]-roots[reals+2*(1-p)]))/2
                    for p,j in enumerate(range(reals,len(roots),2))]
    else:
        numbers += [sum(js[reals:])+sum(mp.sign(roots[reals]-x) for x in roots[:reals])/2]
    if unwrapped:
        return numbers
    return numbers[:reals]+[x-n*mp.nint(x/n) for x in numbers[reals:]]


def quantum_number_error(n, family, holes, branches, roots):
    lengths,labels,_ = configuration(n,family,holes,branches)
    numbers = recovered_quantum_numbers(n,family,holes,branches,roots)
    reals = lengths.count(1)
    errors = [abs(actual-expected) for actual,expected in zip(numbers[:reals],labels[:reals])]
    for actual,expected in zip(numbers[reals:],labels[reals:]):
        delta = actual-expected
        errors.append(abs(delta-n*mp.nint(delta/n)))
    return max(errors)


def solve(n, family, holes, branches):
    if mp.mp.dps < 60:
        raise ValueError('reference requires at least 60 decimal digits')
    lengths, labels, momentum = configuration(n,family,holes,branches)
    guess = [a*mp.tan(mp.pi*i/n)/2 for a,i in zip(lengths,labels)]
    centres = list(mp.findroot(lambda *x: takahashi(n,lengths,labels,x), guess,
                              J=lambda *x: takahashi(n,lengths,labels,x,True),
                              tol=mp.mpf(10)**(-mp.mp.dps+15), maxsteps=100))
    real_count = lengths.count(1)
    reals = centres[:real_count]
    nonself = nonself_conjugate_sign(n,family,labels)
    if nonself:
        # Deliberately low-accuracy seed, not a stored solution or weight.
        # Resolve both conjugate pairs anew at the caller's precision.
        variables = [nonself*mp.mpf('.18'), nonself*mp.mpf('.445'), mp.mpf('.02'),
                     nonself*mp.mpf('.492'), mp.mpf('.96')]
        real_count = 1
    elif family == '22':
        variables = list(reals)
        for i,v in enumerate(centres[real_count:]):
            other = centres[real_count+1-i]
            q = bethe_residual(n,v+.5j,reals+[other+.5j,other-.5j])+1
            variables += [v, mp.mpf('.5')+mp.re(1/(q-1))]
    else:
        v = centres[-1]
        q = bethe_residual(n,v+1j,reals+[v-1j])+1
        displacement = 2j/(q-1)
        variables = reals+[v, v+mp.re(displacement), 1+mp.im(displacement)]
        real_count += 1
    solution = list(mp.findroot(lambda *x: resolved_equations(n,x,real_count), variables,
                               J=lambda *x: resolved_equations(n,x,real_count,True),
                               tol=mp.mpf(10)**(-mp.mp.dps+15), maxsteps=100))
    roots = unpack(solution,real_count)
    if any(solution[j] <= 0 for j in range(real_count+1,len(solution),2)):
        raise ValueError('collapsed or reversed pair')
    if min(abs(x-y) for i,x in enumerate(roots) for y in roots[:i]) < mp.mpf('1e-25'):
        raise ValueError('coincident or unresolved roots')
    residual = max_residual(n,roots)
    if residual > mp.mpf('1e-30'):
        raise ValueError('full rational Bethe residual failed')
    phase = mp.fprod((x+.5j)/(x-.5j) for x in roots)
    if abs(phase-mp.exp(2j*mp.pi*momentum/n)) > mp.mpf('1e-28'):
        raise ValueError('solution left expected momentum sector')
    label_error = quantum_number_error(n,family,holes,branches,roots)
    if label_error > mp.mpf('1e-28'):
        raise ValueError('solution left requested quantum numbers')
    return dict(roots=roots, variables=solution, real_count=real_count,
                centres=centres, momentum_index=momentum, energy=energy(n,roots), residual=residual,
                label_error=label_error, bethe_quantum_numbers=bethe_quantum_numbers(n,roots),
                recovered_labels=recovered_quantum_numbers(n,family,holes,branches,roots),
                recovered_labels_unwrapped=recovered_quantum_numbers(n,family,holes,branches,roots,unwrapped=True),
                root_geometry='nonself-conjugate-3-plus-1' if nonself else 'individual-self-conjugate',
                label_mapping='Deguchi-Giri 2015 Table I' if nonself else 'Hagemans-Caux 2007 Eqs. 25/31')


def fixture_requests():
    cases = []
    for n in (12, 16, 20):
        m = n//4
        cases += [(n, '22', [-m, -m+1, -1, 0, 1, m], ['-.5', '1.5']),
                  (n, '3', [str(-m-.5), str(-m+.5), str(-m+1.5), '-.5',
                            str(m-.5), str(m+.5)], [1])]
    cases += [(12, '3', ['-3.5', '-2.5', '-1.5', '-.5', '.5', '3.5'], [1]),
              (12, '3', ['3.5', '2.5', '1.5', '.5', '-.5', '-3.5'], [-1])]
    return cases


def verify(digits):
    from pathlib import Path
    import re
    if digits < 80:
        raise ValueError('fixture verification requires at least 80 digits')
    fixture = Path(__file__).resolve().parents[1]/'tests/xxx_six_spinon_fixture.hpp'
    chunks = re.findall(r'(?:"[-+0-9.eE]+"\s*)+', fixture.read_text())
    strings = [''.join(re.findall(r'"([-+0-9.eE]+)"', chunk)) for chunk in chunks]
    with mp.workdps(digits):
        expected = []
        for n, family, holes, branches in fixture_requests():
            state = solve(n, family, holes, branches)
            expected += [state['energy'], raising_weight(n, ground(n), state['roots'])]
            expected += [value for z in state['roots'] for value in (mp.re(z), mp.im(z))]
        if len(strings) != len(expected):
            raise RuntimeError('unexpected C++ fixture structure')
        accuracy = min(90, digits-25)
        for text, value in zip(strings, expected):
            tolerance = mp.mpf(10)**(-accuracy) * (abs(value) if value else 1)
            if abs(mp.mpf(text)-value) > tolerance:
                raise RuntimeError('C++ fixture differs from independent reference')
    print(f'Eight selected six-spinon root/energy/weight fixtures agree to {accuracy} digits')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sites', type=int, nargs='?')
    parser.add_argument('family', choices=['22', '3'], nargs='?')
    parser.add_argument('--holes', nargs=6)
    parser.add_argument('--branches', nargs='+')
    parser.add_argument('--digits', type=int, default=110)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify(args.digits)
        return
    if args.sites is None or args.family is None or args.holes is None or args.branches is None:
        parser.error('supply sites, family, six holes and string branches, or use --verify')
    with mp.workdps(args.digits):
        state = solve(args.sites, args.family, args.holes, args.branches)
        mu = ground(args.sites)
        weight = raising_weight(args.sites, mu, state['roots'])
        digits = args.digits-20
        print(json.dumps(dict(sites=args.sites, family=args.family, holes=args.holes, branches=args.branches,
            digits=args.digits, momentum_index=state['momentum_index'],
            energy=mp.nstr(state['energy'], digits), gap=mp.nstr(state['energy']-energy(args.sites,mu), digits),
            raising_weight=mp.nstr(weight, digits), residual=mp.nstr(state['residual'], 8),
            label_error=mp.nstr(state['label_error'], 8), root_geometry=state['root_geometry'],
            label_mapping=state['label_mapping'], recovered_labels=[mp.nstr(x,digits) for x in state['recovered_labels']],
            roots=[[mp.nstr(mp.re(z),digits),mp.nstr(mp.im(z),digits)] for z in state['roots']]), indent=2))


if __name__ == '__main__':
    main()
