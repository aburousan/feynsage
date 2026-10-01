# feynsage in detail

The front page (`README.md`) is short. This page keeps the full description of every module, every check, the speed numbers and the limits.

## What it does

| Module | Content |
|---|---|
| `graph.FeynmanGraph` | A diagram as a multigraph with external momenta. Spanning trees and 2-forests are enumerated from their definitions and give the Symanzik polynomials U and F. U is also computed from the Kirchhoff matrix-tree theorem as an independent check. `family()` routes the momenta through a cycle basis. |
| `family.IntegralFamily` | Propagators, scalar products, U and F by the matrix method (U = det M), Lee's zero-sector criterion and all L(L+E) IBP identities for any index vector. |
| `laporta.Reducer` | Laporta reduction with exact rational functions in d and the invariants, user-given symmetries and the list of master integrals. |
| `ff.reduce_ff` | The same reduction done the way FIRE and Kira do it: the IBP identities are written once with symbolic indices, the system is solved numerically modulo large primes, only the equations the targets need are kept after the first probe and the coefficients are rebuilt by Thiele interpolation in d (and one invariant) and rational reconstruction. |
| `pv` | One-loop tensor integrals in the Package-X / LoopTools normalisation. `loop(numerator, [q, m], ...)` takes strings (l^mu, l.p, l^2, g(mu,nu)) and reduces to A0, B0, C0, D0 by projection on a symmetric basis, with the Gram matrix inverted exactly (adjugate over the polynomial ring) and all arithmetic in the field Q(d, invariants, masses). Rational terms come from d = 4 - 2 eps. Light-like two-point tensors (p^2 = 0) are done by Feynman parameters without a Gram matrix. `PVB`, `PVC`, `PVD` as Package-X's. A0 and B0 in closed form with the +i0 branch (`LogM`, `DiscB`). IR-finite C0 and D0 stay as symbols whose `.n()` gives more than 30 digits, and `explicit()` writes them in dilogarithms (`scalar.c0_closed`, `scalar.d0_closed`). If the kinematics make the Gram determinant vanish (all invariants zero, g-2 at q^2 = 0), the tensor coefficients come straight from Feynman parameters, exactly. `euclidean=True` for Euclidean integrals. |
| `scalar` | Closed forms of C0 and D0 for real kinematics: logarithms and dilogarithms of exact algebraic numbers, every branch fixed by the +i0. C0 from its one-dimensional representation (Denner's twelve dilogarithms for lambda > 0, explicit branch bookkeeping for lambda < 0). D0 from Denner (1993) eq. (4.43) when no mass vanishes and from Denner and Dittmaier (2010) eqs. (3.76) to (3.84) for one to four massless lines. `c0_value`, `d0_value` give the numbers. |
| `ir` | IR-divergent (soft or collinear) C0 and D0 in dimensional regularisation: the six triangles and sixteen boxes of Ellis and Zanderighi (2008), with every relabelling found automatically, written in Package-X's normalisation. `C0(...)` and `D0(...)` return c_-2/eps^2 + c_-1/eps + c_0 with mu for these, and `loop()`, `PVC`, `PVD` carry the poles through the reduction, keeping the d-dependence of the coefficients to second order. |
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
- `test_pv.sage`: A0, B0, B1, B00, B11, B0', C0, C1, C2, C00, C11, C12 and D0 against Package-X 2.1.1 at 190 kinematic points (below, between and above thresholds, light-like momenta, massless lines, vanishing Gram determinants): A, B and the C tensors to 1e-15; C0 and D0 agree with every digit of Package-X's double-precision values. The reference values are in `tests/data`.
- `test_ir.sage`: all six IR-divergent triangles and sixteen boxes, each at a random relabelling as well. The 1/eps^2, 1/eps and eps^0 parts agree with Package-X to double precision at 118 points. Boxes 15 and 16 (Package-X does not evaluate their finite parts) and the threshold points are checked against a photon-mass regularised D0.
- `test_d0.sage`: D0 against Package-X at 146 points (zero to four massless lines) to 1e-9, against itself under all 24 relabellings of the lines to 1e-34, and against direct integration at the two points where Package-X itself is only good to 1e-8.
- `test_ir_tensor.sage`: tensor coefficients with IR poles (PVC, PVD) and at vanishing Gram determinants against Package-X. One Package-X value is wrong: its numerical C0(-4,-3,-1; 0,0,1) = +0.671, but all invariants are spacelike, so C0 must be negative. Direct integration gives -0.585977, as feynsage does, and so do Mathematica's NIntegrate and FeynCalc's own Feynman parametrisation integrated numerically. Package-X is negative and smooth at nearby points and gives a third value (-0.721) for the same integral with the lines relabelled. LoopTools 2.16 fails at the same point: its FF version gives the same +0.671 and -0.721, and its Denner version gives complex numbers that change with the labelling (`tests/looptools`). The test checks this.
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

In `pv`: the closed forms of C0 and D0 need numerical kinematics. With symbols, an IR-finite
C0 or D0 stays a symbol. For IR-divergent ones whose formula needs the roots K or gamma of Ellis
and Zanderighi (triangle 6, boxes 11 to 16), the coefficients become the named functions
`C0IR(k, ...)`, `D0IR(k, ...)`, which evaluate once numbers are put in. Points where no
relabelling gives a regular closed form fall back to contour integration. At these points the
Cayley determinant vanishes, a leading Landau singularity, where the box itself is singular.
The Feynman-parameter route for a vanishing Gram determinant covers one group of lines with
identical kinematics (any masses), or two groups with equal masses inside each group. Other
cases raise an error; move off that point.


## Install, in detail

`./install.sh` works on macOS (Apple silicon and Intel), on Linux (x86_64 and arm64) and on
Windows through WSL. It finds SageMath (the macOS app, a conda environment or a `sage` on the
PATH) and installs it if it is missing (Homebrew cask on macOS, a conda-forge environment
`sage` on Linux). It installs FORM if it is missing (Homebrew on macOS, the static release
binary from github.com/form-dev/form into `~/.local/bin` on Linux) and then installs feynsage
into Sage's Python. On the macOS app this is a user install, so the signed app is never
modified. A C compiler is optional: with one, the finite-field reducer compiles its kernel on
first use; without one it runs in pure Python with the same answers.

Tested on macOS (Apple silicon) with SageMath 10.9 and on CentOS 7 (glibc 2.17, GCC 4.8) with
SageMath 10.7 from conda-forge.
