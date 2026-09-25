# Periodic Hubbard filling, spin sectors, and attractive U

[Hubbard equations and numerical method](hubbard.md) · [Overview](../README.md)

Use `--particles N` and `--sz VALUE` to select the physical sector. The two
spin populations are $`N_{\mathrm{up}} =N /2+S^z`$ and $`N_{\mathrm{down}} =N /2-S^z`$; each must be an integer
between zero and L. Half-integers can be written as `1/2` or `0.5`.
The defaults remain N=L and Sz=0, and L must be even.
This guide describes `bethe-hubbard-pbc`. The [free-end solver](hubbard-open.md)
uses the same energy mappings without the periodic shell-parity restrictions;
it also accepts odd L and defaults to Sz=1/2 when N is odd.

## Supported ground states

The numerical core solves a repulsive sector with N<=L and N_down<=N_up.
Exact symmetries reduce other physical sectors to this core.

| Physical request | Current support |
| --- | --- |
| U=0 | Every valid N and Sz; independent free-fermion bands |
| Only one spin species | Every valid N and either sign of U; no double occupation |
| U>0, N=L | Every valid Sz |
| U>0, N<L | Odd N_up and odd N_down, plus the single-species case |
| U>0, N>L | Full particle-hole image of a supported sector below half filling |
| U<0, Sz=0 | Every even N from 0 to 2L, via a repulsive half-filled spin sector |
| U<0, Sz!=0 | Supported when the mapped repulsive sector belongs to the families above |

For example:

```sh
# Repulsive doping: three electrons of each spin.
build/bethe-hubbard-pbc 16 --u 4 --particles 6
# Electron doping: maps to six particles below half filling.
build/bethe-hubbard-pbc 16 --u 4 --particles 26
# Balanced attraction: four pairs, mapped to polarized half filling.
build/bethe-hubbard-pbc 16 --u -4 --particles 8 --roots
# A half-filled polarized repulsive sector.
build/bethe-hubbard-pbc 16 --u 4 --sz -2
# Odd N is allowed on the unrestricted free-fermion path.
build/bethe-hubbard-pbc 16 --u 0 --particles 7 --sz 1/2
```

This is deliberately not yet an arbitrary-filling, arbitrary-spin ground-state
solver. In particular, `16 --u 4 --particles 8 --sz 0` is rejected: both spin
populations are even. A periodic ring's shell parity can change which spin
branch is lowest. Continuing a single compact quantum-number sea is not a
general minimization over those branches or their spin descendants. Supporting
these remaining sectors needs branch selection, not just more Newton variables.
No complex-string or general excited-state solver is implied here.

## Why attraction can use a repulsive solver

On an even bipartite ring, the down-spin particle-hole transformation

```math
c_{j,\downarrow}\longrightarrow(-1)^j c^\dagger_{j,\downarrow}.
```

preserves hopping and exchanges the sign of U. This is the Shiba mapping;
see [rylands-2022](../CITATIONS.md#rylands-2022), Sec. II. For our **unshifted**
Hamiltonian, writing g=|U|,

```math
\begin{aligned}
(N_\uparrow,N_\downarrow)&\longrightarrow(N_\uparrow,L-N_\downarrow),\\{}
E_{-g}(N_\uparrow,N_\downarrow)&=E_{+g}(N_\uparrow,L-N_\downarrow)-gN_\uparrow.
\end{aligned}
```

Thus balanced attraction maps to N'=L with `Sz'=(N-L)/2`. The attractive
energy can be found with real repulsive roots even though the physical
attractive Bethe roots need not be real. This is an exact unitary mapping,
not a string approximation. We return the auxiliary roots explicitly; we
do not reconstruct the physical attractive roots or wavefunction.

Above half filling, a particle-hole transformation on **both** species gives

```math
\begin{aligned}
(N_\uparrow,N_\downarrow)&\longrightarrow(L-N_\uparrow,L-N_\downarrow),\\{}
E_U(N_\uparrow,N_\downarrow)&=E_U(L-N_\uparrow,L-N_\downarrow)+U(N-L).
\end{aligned}
```

Finally, spin reversal exchanges N_up and N_down without changing energy.
The implementation applies Shiba first when U<0, then full particle-hole if
needed, then spin reversal. These operations preserve the minimum within the
mapped sector. Odd rings are excluded because the staggered transformation
would change the periodic boundary condition.

## Momentum and root metadata

Momentum also transforms. In the convention `P=sum(k) mod 2*pi`, for even L,
the down-spin Shiba step adds $`\pi \,(N_{\mathrm{down}} -1)`$ when converting auxiliary
momentum back to the physical sector. The full particle-hole step adds
$`\pi \,N`$. Counts in these expressions belong to the sector **before that step**;
spin reversal adds nothing. These offsets follow by replacing each occupied
down-spin momentum by its complementary hole at $`\pi -k`$; a filled one-spin
band has momentum pi. Tests check the resulting physical momenta directly
against translation in Fock space, including degenerate ground eigenspaces.

In `State<Real>`:

- `particles`, `down_spins`, `interaction`, `energy`, `momentum_index`, and
  `momentum` describe the requested **physical** sector.
- `root_particles`, `root_down_spins`, and `root_interaction` describe
  `charge_momenta`, `spin_rapidities`, `quantum_numbers`, and all residuals.
- `energy_offset` and `momentum_offset` obey `E=E_roots+energy_offset` and
  `momentum_index=(root_index+momentum_offset)%L`.
- `shiba_transformed`, `particle_hole_transformed`, and `spin_reversed` record
  the mapping; `auxiliary_roots()` indicates whether any was used.

The CLI prints this distinction whenever a mapping was applied, including
with `--roots` and for unconverged estimates. Free-fermion occupation lists
have no Bethe labels or spin rapidities. For open shells, their deterministic
representative fills `q=-(count-1)/2,...` with integer division toward zero;
another degenerate representative may have a different momentum.

Both energy offsets and all solves use the selected scalar precision.
No Uni20 linear-algebra extension is needed: the same real dense solve serves
these sectors. Residuals are equation diagnostics, not energy-error bounds;
in particular, adding an O(|U|N) attractive offset can conceal small binding
or excitation scales at finite precision.
