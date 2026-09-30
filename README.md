# feynsage

Feynman-integral tools written natively for SageMath. It covers the parts of the lecture
notes "Feynman integrals at one loop" that a package like Package-X or LiteRed does,
with Sage's graph theory and exact linear algebra underneath and FORM for Dirac
algebra.

## Install

```bash
git clone https://github.com/aburousan/feynsage.git
cd feynsage
./install.sh          # add --test to run the test suite afterwards
```

The script works on macOS (Apple silicon and Intel), on Linux (x86_64 and arm64) and on
Windows through WSL. It finds SageMath (the macOS app, a conda environment or a `sage` on the
PATH) and installs it if it is missing (Homebrew cask on macOS, a conda-forge environment
`sage` on Linux). It installs FORM if it is missing (Homebrew on macOS, the static release
binary from github.com/form-dev/form into `~/.local/bin` on Linux) and then installs feynsage
into Sage's Python. On the macOS app this is a user install, so the signed app is never
modified. A C compiler is optional: with one, the finite-field reducer compiles its kernel on
first use; without one it runs in pure Python with the same answers.

Tested on macOS (Apple silicon) with SageMath 10.9 and on CentOS 7 (glibc 2.17, GCC 4.8) with SageMath 10.7
from conda-forge.

## Walkthrough notebook

`examples/feynsage_walkthrough.ipynb` goes through the package step by step with a diagram
for every step: momenta, spanning trees and 2-forests, U and F three ways, zero sectors
(drawn by shrinking lines), IBP identities, Laporta and finite-field reduction, one-loop
formulas, FORM traces, plots and at the end a real prediction: the lifetime of the neutral
pion from the quark triangle, 8.35e-17 s against the measured (8.43 +- 0.13)e-17 s.

## What it does

| Module | Content |
|---|---|
| `graph.FeynmanGraph` | A diagram as a multigraph with external momenta. Spanning trees and 2-forests are enumerated from their definitions and give the Symanzik polynomials U and F. U is also computed from the Kirchhoff matrix-tree theorem as an independent check. `family()` routes the momenta through a cycle basis. |
| `family.IntegralFamily` | Propagators, scalar products, U and F by the matrix method (U = det M), Lee's zero-sector criterion and all L(L+E) IBP identities for any index vector. |
| `laporta.Reducer` | Laporta reduction with exact rational functions in d and the invariants, user-given symmetries and the list of master integrals. |
| `ff.reduce_ff` | The same reduction done the way FIRE and Kira do it: the IBP identities are written once with symbolic indices, the system is solved numerically modulo large primes, only the equations the targets need are kept after the first probe and the coefficients are rebuilt by Thiele interpolation in d (and one invariant) and rational reconstruction. |
| `oneloop` | Closed results: tadpole, massless bubble with any powers, on-shell triangle with any powers, equal-mass and one-mass bubbles and the epsilon expansion (Gamma poles are made explicit before expanding). |
| `form` | Dirac traces with FORM: `trace(indices)` in D dimensions and `dirac_trace(factors, vectors)` for products with slashed momenta, masses and gamma_5 (for example the pion triangle, `4*m*eps(k1, k2, mu, nu)`). FORM is found on the PATH, in `~/.local/bin`, `~/bin`, Homebrew, or through `FEYNSAGE_FORM`. |
| `plotting` | Figures in the style of the notes: `set_theme()` (also for Sage's own `plot`), `curves`, `mb_plane`, `FeynmanGraph.plot()` (vertices on a circle with fewest crossings, or on a line for chains; labels x_i matching U and F; thick blue massive lines; pink dots for raised powers), `draw_panels` (spanning trees and 2-forests, removed lines dashed) and `draw_sectors` (sectors, with the lines of index 0 shrunk to points by `FeynmanGraph.contract`). `sage examples/plots.sage` writes examples. |

## Checks

`tests/` reproduces results that were verified independently in the notes:

- `test_bubble.sage`: the equal-mass bubble reductions agree with Kira 3.1.
- `test_kite.sage`: the two-loop propagator family: masters, F(1,1,1,1,1) and F(2,2,1,2,2) agree with LiteRed.
- `test_graphs.sage`: U and F of the kite, the massless box, the two-loop phi^4 diagram and the two-mass bubble agree three ways (tree rules, Kirchhoff, matrix method).
- `test_ff.sage`: the finite-field reducer gives the same kite and bubble results (LiteRed, Kira).
- `test_oneloop_form.sage`: Lee's criterion on massless and massive tadpoles, the massless triangle, the two-loop phi^4 diagram of the notes and a FORM trace.

Run them with `sage tests/<name>.sage` from this folder.

## Example

```python
from feynsage import *
kin = Kinematics(['p'], ['s'], {('p','p'): 's'}, euclidean=False)
fam = IntegralFamily('bub', ['l'], kin, [(mom(l=1), 1), (mom(l=1, p=-1), 1)])
red = Reducer(fam, symmetries=[(1, 0)]).run(rmax=3)
red.masters()          # [(0, 1), (1, 1)]
red.reduce((2, 1))     # {(1,1): (-d+3)/(s-4), (0,1): (d-2)/(2*s-8)}

# the same with finite fields, much faster on bigger systems
red = Reducer(fam, symmetries=[(1, 0)])
reduce_ff(red, [(2, 1), (3, 1)], rmax=3)

g = FeynmanGraph([('L','B',0), ('L','T',0), ('B','R',0), ('T','R',0), ('T','B',0)],
                 {'L': mom(p=1), 'R': mom(p=-1)},
                 Kinematics(['p'], ['pp'], {('p','p'): 'pp'}, euclidean=True))
g.U(), g.F()           # the kite, read off its trees and 2-forests
```

## Practical points

- Give `run(rmax, smax)` numerator seeds (`smax >= 1`) whenever the family has more than one loop: IBP identities produce integrals with numerators and they can only be reduced if seeds with numerators are included.
- Symmetries can be given as generators; the reducer closes them into a group.
- If the problem has a single scale, set it to 1 in the kinematics rules (like Kira's `symbol_to_replace_by_one`). The arithmetic is then univariate in d and runs in seconds; the scale comes back by dimensional analysis. The two-loop test goes from more than ten minutes to about 7 seconds this way.

## Speed

Two-loop kite family (massless, q^2 = 1), coefficients checked to be identical in every case:

| Seeds | Equations | Exact `Reducer` | `reduce_ff` | Kira 3.1 |
|---|---|---|---|---|
| rmax 4, smax 2 | 9 368 | 7.3 s | 0.2 s | |
| rmax 9, smax 4 (target kite[3,3,2,3,3]) | 157 843 | too slow | 2.7 s | 6.7 s |

feynsage ran on the laptop, Kira on a server with one thread, so the last column is only a rough guide.

Where the time goes in the rmax 9 run: building the equations with numpy (all seeds of one IBP
template at once) 1.0 s; the first probe, a sparse elimination of all 157 843 equations in C
(`csrc/elim.c`, compiled with the system `cc` on first use and loaded with ctypes) 0.9 s; it
keeps only the 895 equations the targets need and every later probe on those takes 8 ms in
Python. Without a C compiler the same steps run in pure Python, about ten times slower.

## Limits

`reduce_ff` handles d and at most one kinematic invariant, so any other scale must be set to a
number. The seeds are generated for whole sectors without Kira's selection of seeds, so very
large systems need more memory than Kira. Neither reducer finds symmetries by itself (they can
be given them). The seed range is chosen by hand: `rmax` counts dots and `smax` numerator
powers (Kira's `r` is the total of the positive powers instead) and a target outside the
range raises an error.

## Licence

MIT, see `LICENSE`.
