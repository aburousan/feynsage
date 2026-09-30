# feynsage

Feynman-integral tools written natively for SageMath. It does what Package-X does at one
loop (tensor reduction to A0, B0, C0, D0 in the same conventions, checked point by point
against it) and what LiteRed, FIRE and Kira do beyond (graph polynomials, IBP reduction
with finite fields), with Sage's graph theory and exact algebra underneath and FORM for
Dirac algebra. It grew out of the lecture notes "Feynman integrals at one loop".

## Quick start

```python
from feynsage import *

# one loop, like Package-X: exact, symbols allowed
loop("l^mu l^nu", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})   # g^{mu nu} B00 + p^mu p^nu B11
B0(s, m, m)                                                           # closed form with the +i0 branch
C0(0, 0, 5, 1, 1, 1).n()                                              # numbers when you want them

# any number of loops: families from strings, reduction with automatic seeds and symmetries
kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
ibp_reduce(kite, ["F(1,1,1,1,1)", "F(2,2,1,2,2)"])

# graphs: U and F from spanning trees and 2-forests and a picture
g = diagram("kite"); g.U(), g.F(); g.plot()

# every function explains its output
loop("l.p", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"}, explain=True)
info(ibp_reduce)
```

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

## Notebooks

- `examples/feynsage_walkthrough.ipynb` goes through the package step by step with a diagram
  for every step: momenta, spanning trees and 2-forests, U and F three ways, zero sectors
  (drawn by shrinking lines), IBP identities, Laporta and finite-field reduction, one-loop
  formulas, FORM traces, plots, a real prediction (the lifetime of the neutral pion from the
  quark triangle, 8.35e-17 s against the measured (8.43 +- 0.13)e-17 s), one-loop tensor
  integrals the Package-X way, the four examples of the lecture notebook and the short front end.
- `examples/peskin_examples.ipynb` does two QED classics from Peskin and Schroeder: the vacuum
  polarisation (exact Ward identity, eq. (7.91), running of alpha: the leptons give
  Delta alpha(M_Z^2) = 0.031419 at one loop) and the anomalous magnetic moment, F2(0) = alpha/2pi
  exactly, with U and F from the vertex graph and the Dirac algebra in FORM.

## What it does

| Module | Content |
|---|---|
| `graph.FeynmanGraph` | A diagram as a multigraph with external momenta. Spanning trees and 2-forests are enumerated from their definitions and give the Symanzik polynomials U and F. U is also computed from the Kirchhoff matrix-tree theorem as an independent check. `family()` routes the momenta through a cycle basis. |
| `family.IntegralFamily` | Propagators, scalar products, U and F by the matrix method (U = det M), Lee's zero-sector criterion and all L(L+E) IBP identities for any index vector. |
| `laporta.Reducer` | Laporta reduction with exact rational functions in d and the invariants, user-given symmetries and the list of master integrals. |
| `ff.reduce_ff` | The same reduction done the way FIRE and Kira do it: the IBP identities are written once with symbolic indices, the system is solved numerically modulo large primes, only the equations the targets need are kept after the first probe and the coefficients are rebuilt by Thiele interpolation in d (and one invariant) and rational reconstruction. |
| `pv` | One-loop tensor integrals in the Package-X / LoopTools normalisation. `loop(numerator, [q, m], ...)` takes strings (l^mu, l.p, l^2, g(mu,nu)) and reduces to A0, B0, C0, D0 by projection on a symmetric basis, with the Gram matrix inverted exactly (adjugate over the polynomial ring) and all arithmetic in the field Q(d, invariants, masses). Rational terms come from d = 4 - 2 eps. Light-like two-point tensors (p^2 = 0) are done by Feynman parameters without a Gram matrix. `PVB`, `PVC` as Package-X's. A0 and B0 in closed form with the +i0 branch (`LogM`, `DiscB`); C0 and D0 numerically for any real kinematics by contour deformation of the Feynman-parameter integral, with a two-grid error check. `euclidean=True` for Euclidean integrals. |
| `easy` | `family(props, kin)`, `graph("A-B, A-B:m", legs, kin)`, `diagram(name)`, `symmetries(fam)` (automorphisms of U + F found with Sage's graph automorphism code: Pak's criterion) and `ibp_reduce(fam, targets)`, which picks the seed range and the reducer. |
| `explain` | `info(f)` and `explain=True` on every public function print what it computes and what its output means. |
| `oneloop` | Closed results: tadpole, massless bubble with any powers, on-shell triangle with any powers, equal-mass and one-mass bubbles and the epsilon expansion (Gamma poles are made explicit before expanding). |
| `form` | Dirac traces with FORM: `trace(indices)` in D dimensions and `dirac_trace(factors, vectors)` for products with slashed momenta, masses and gamma_5 (for example the pion triangle, `4*m*eps(k1, k2, mu, nu)`). FORM is found on the PATH, in `~/.local/bin`, `~/bin`, Homebrew, or through `FEYNSAGE_FORM`. |
| `plotting` | Figures in the style of the notes: `quick_plot(expr, (s, a, b))` (one line, real and imaginary parts, finite part of expressions with 1/eps, parallel evaluation), `set_theme()` (also for Sage's own `plot`), `curves`, `mb_plane`, `FeynmanGraph.plot()` (vertices on a circle with fewest crossings, or on a line for chains; labels x_i matching U and F; thick blue massive lines; pink dots for raised powers), `draw_panels` (spanning trees and 2-forests, removed lines dashed) and `draw_sectors` (sectors, with the lines of index 0 shrunk to points by `FeynmanGraph.contract`). `sage examples/plots.sage` writes examples. |

## Checks

`tests/` reproduces results that were verified independently in the notes:

- `test_bubble.sage`: the equal-mass bubble reductions agree with Kira 3.1.
- `test_kite.sage`: the two-loop propagator family: masters, F(1,1,1,1,1) and F(2,2,1,2,2) agree with LiteRed.
- `test_graphs.sage`: U and F of the kite, the massless box, the two-loop phi^4 diagram and the two-mass bubble agree three ways (tree rules, Kirchhoff, matrix method).
- `test_ff.sage`: the finite-field reducer gives the same kite and bubble results (LiteRed, Kira).
- `test_oneloop_form.sage`: Lee's criterion on massless and massive tadpoles, the massless triangle, the two-loop phi^4 diagram of the notes and a FORM trace.
- `test_pv.sage`: A0, B0, B1, B00, B11, B0', C0, C1, C2, C00, C11, C12 and D0 against Package-X 2.1.1 at 190 kinematic points (below, between and above thresholds, light-like momenta, massless lines): A and B to 1e-15, C to 1e-11, D0 to 1e-6. The reference values are in `tests/data`.
- `test_peskin.sage`: vacuum polarisation (Ward identity exactly, Peskin (7.91) to 1e-15, one-loop leptonic running of alpha) and F2(0) = alpha/(2 pi) exactly.
- `test_professor_uf.sage`: U and F of the four examples of the lecture notebook (triangle, massive bubble, box, two-loop sunset), both from the propagators and from the graph.

Run them with `sage tests/<name>.sage` from this folder.

## Example with the lower-level classes

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
- `ibp_reduce` finds the symmetries itself; with `Reducer` directly they can be given as generators and the reducer closes them into a group.
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
`nproc > 1` (in `reduce_ff` and `ibp_reduce`) evaluates the sample points in parallel worker
processes and gives the same result. On hercules (an older 192-core server) the rmax 9 kite goes
from 11.6 s with one process to 8.1 s with four and 7.9 s with eight: the build and the first
probe (about 5 s there) are still serial, so the gain grows with the size of the trimmed system
(`tests/bench_parallel.sage`).

One-loop tensor reduction with everything symbolic (all masses and invariants symbols, on the laptop):
rank-2 bubble 1.1 s, rank-1 triangle 1.6 s, C00 0.06 s, rank-3 triangle 4.7 s, rank-2 box 1.9 s.
With numbers it is faster still; `quick_plot(..., nproc=4)` evaluates plot points in parallel.

## Limits

`reduce_ff` handles d and at most one kinematic invariant, so any other scale must be set to a
number. The seeds are generated for whole sectors without Kira's selection of seeds, so very
large systems need more memory than Kira. With the low-level functions the seed range is chosen
by hand: `rmax` counts dots and `smax` numerator powers (Kira's `r` is the total of the
positive powers instead) and a target outside the range raises an error.

In `pv`: C0 and D0 are exact symbols in every formula but their values are numerical (Package-X
has closed forms in dilogarithms). IR-divergent C0 and D0 (massless internal lines between
on-shell legs) are not separated into poles; the eps expansion assumes they are IR finite.
Kinematics with a vanishing Gram determinant for three or four points (for example all
invariants zero, or collinear momenta) raise an error; move off that point or keep symbols.

## Licence

MIT, see `LICENSE`.
