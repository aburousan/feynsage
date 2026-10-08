# feynsage

Feynman integrals in SageMath, written from scratch.

- **One loop, like Package-X.** Tensor integrals reduced to A0, B0, C0, D0, exactly, in the same
  conventions. C0 and D0 in closed form (logarithms and dilogarithms). Soft and collinear poles
  written out. Works also when the Gram determinant is zero. Raised propagator powers, Taylor
  series, discontinuities and double spectral functions, pentagons, fermion lines with the Dirac
  equation and the Gordon identities, form-factor projectors, transverse and longitudinal parts.
- **Many loops, like LiteRed, FIRE and Kira.** Integral families from strings, IBP reduction
  (exact, or with finite fields), symmetries found by itself, master integrals.
- **Graphs.** The Symanzik polynomials U and F from spanning trees and 2-forests, checked against
  two other methods. Pictures of diagrams, trees and sectors.
- **Dirac algebra with FORM**, gamma_5 included. Traces, spin sums, polarization sums and colour factors
  are written in textbook notation, `dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))`, and FORM
  does the work without any FORM code ([docs/dirac_algebra.md](docs/dirac_algebra.md)).
- **Diagrams from a model, like FeynArts.** `topologies`, `insert_fields` with the Standard Model of
  FeynArts (`SM()`, also QED), pictures of the diagrams, and tree amplitudes squared with FORM, with
  Breit-Wigner widths if wanted.
- **The methods of the lecture note**: Laurent series at Gamma poles, Mellin-Barnes residues, differential
  equations for master integrals, Goncharov polylogarithms, sector decomposition, expansion by regions,
  counting master integrals.
- Every function can explain its own output: `explain=True` or `info(f)`.

I wrote feynsage in 2024 for my MSc project. After the talk "An Introduction to Feynman Integrals" by
Prof. B. Ananthanarayan (NISER, 2026) it was refined, checked against the lecture notes of that
talk and published. His introductory lectures on the same topic at IMSc, Chennai, are on YouTube:
[Introductory lectures on Feynman Integrals Part1](https://youtu.be/tvCBzU9GWeI).
A step-by-step tutorial with real computations is in [`docs/tutorial/feynsage_tutorial.pdf`](docs/tutorial/feynsage_tutorial.pdf).

## Who it is for

feynsage is meant for small and medium problems and for studying: every step is exact, can be
checked and explains itself. On such problems it runs in the same range as the established tools
and sometimes a little faster. For very large reductions (millions of equations, several scales)
Kira and FIRE remain the better choice: they are built for that scale and are much more economical
with memory. feynsage uses all cores for the parts that parallelise (sample points of the
finite-field reduction, many one-loop numbers at once, plots). I may continue improving feynsage in the future.

## Install

```bash
git clone https://github.com/aburousan/feynsage.git
cd feynsage
./install.sh          # add --test to run all the tests afterwards
```

macOS, Linux and Windows (WSL). The script finds or installs SageMath and FORM.

## Try it

```python
from feynsage import *

# one loop: exact, symbols allowed
loop("l^mu l^nu", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})   # g^{mu nu} B00 + p^mu p^nu B11
C0(0, 0, 5, 1, 1, 1).n(digits=30)                                   # 30 digits, all correct
explicit(D0(1, 2, 3, 4, -5, -6, 1, 1, 1, 1))                        # D0 in dilogarithms
C0(0, 5, 0, 0, 0, 0)                                                # IR divergent: 1/eps^2, 1/eps shown

# two loops: reduce the kite to master integrals
kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
ibp_reduce(kite, ["F(2,2,1,2,2)"])                       # finite fields
ibp_reduce(kite, ["F(2,2,1,2,2)"], method="trimmed")     # exact arithmetic, just as fast

# graph polynomials and a picture
g = diagram("kite"); g.U(), g.F(); g.plot()
```

## Function reference

[`docs/REFERENCE.md`](docs/REFERENCE.md) lists every public function with its arguments and what it
returns. `info(f)` prints the same for one function inside Sage.

## Notebooks

The notebooks in `examples/` are written for someone who knows only a little Python (`print`, functions and
`for` loops). Run them in order with the SageMath kernel; every output in them is real.

- [`examples/00_start_here.ipynb`](examples/00_start_here.ipynb): Python and Sage from zero: variables, lists,
  loops, functions, symbols, simplifying, calculus and plots.
- [`examples/01_qft_tree_level.ipynb`](examples/01_qft_tree_level.ipynb): a first course in QFT. Dirac traces and
  spin sums by hand, then `process(...)` for QED, QCD and the electroweak theory: |M|^2, helicity amplitudes,
  cross sections, decay widths, e+e- -> W+W- and its gauge cancellation, 2 -> 3. Checked against Peskin and
  Schroeder (5.11)-(5.13), (5.21), (5.87), (5.105)-(5.106), (21.107) and the PDG top width.
- [`examples/02_qft_one_loop.ipynb`](examples/02_qft_one_loop.ipynb): one loop. A0, B0, C0, D0 and tensor
  reduction, the photon and electron self-energies (Peskin (7.90), (7.91), (10.41)), the running of alpha,
  F2(0) = alpha/2pi in two gauges, the QCD counterterms and the beta function (Peskin (16.73)-(16.85)), H -> gg
  and H -> gamma gamma, and the one-loop correction to a cross section with `virtual()`.
- [`examples/03_lecture_note_feynman_integrals.ipynb`](examples/03_lecture_note_feynman_integrals.ipynb): the
  lecture note *Feynman integrals at one loop* of Prof. B. Ananthanarayan (NISER, 2026), section by section with
  its feynsage boxes: the tadpole, dimensional regularisation, IBP, the bubble and its differential equation,
  Feynman parameters and graph polynomials, Mellin-Barnes, Laporta, master integrals, canonical differential
  equations, sector decomposition, expansion by regions, Kira and the decay of the neutral pion.
- [`examples/04_feynsage_toolbox.ipynb`](examples/04_feynsage_toolbox.ipynb): the tools one by one: graphs,
  U and F, families, IBP reduction (exact and finite fields), master counting, differential equations, sector
  decomposition, Mellin-Barnes, the Package-X layer, FORM, plots, parallel work and the built-in help.
- [`docs/tutorial/feynsage_tutorial.pdf`](docs/tutorial/feynsage_tutorial.pdf): the step-by-step tutorial in
  the style of the lecture note, with the code, its real output and an explanation of each function used.
- `examples/chluba/`: Compton and double Compton scattering from J. Chluba's thesis (CMB spectral
  distortions) with feynsage and with FeynCalc, and the cancellation of the infrared divergence in
  dimensional regularisation.

## Is it right?

The tested results are checked against another program, a known answer or a direct integral:

| Part | Checked against | Agreement |
|---|---|---|
| A0, B0, C and D tensors | Package-X 2.1.1, 190 points | every digit of its values |
| IR-divergent C0, D0 | Package-X, 118 points; a photon mass where Package-X has no number | double precision |
| D0, also with massless lines | Package-X (146 points), LoopTools, direct integration | 1e-9 (Package-X's own limit), 1e-15 against the others |
| IBP reduction | Kira 3.1 (80 runs: kite, vertex, sunset, double box, see [`benchmarks/results.md`](benchmarks/results.md)), LiteRed | identical, coefficient by coefficient |
| U and F | sir's notebook (uf-new-short.nb); FeynCalc | identical |
| QED | Peskin and Schroeder | Ward identity and F2(0) = alpha/2pi exactly |
| Series, discontinuities, raised powers, fermion lines, projectors | Package-X 2.1.1 (`tests/test_px_parity.sage`, `tests/test_dirac_px.sage`); finite differences where Package-X gives Indeterminate | every digit printed |
| Pentagons | LoopTools 2.16 (E0), direct Feynman-parameter integration | 1e-14 |
| Tree-level \|M\|^2 (`process`) | FeynArts + FeynCalc 10.2: 38 SM processes with every mass kept, among them e+e- -> W+W-, ZZ, ZH, t t~, W+W- -> W+W-, gamma gamma -> W+W-, g g -> t t~, u g -> d W+ and decays (`tests/test_feyncalc.sage`); 2 -> 3 and 1 -> 3 processes with masses (e+e- -> mu+mu- gamma, e+e- -> u u~ g, e+e- -> Z H gamma, mu- -> e- nu_e~ nu_mu, t -> b e+ nu_e, H -> e- e+ gamma) at random phase-space points, epsilon tensors included (`tests/test_feyncalc_3body.sage`) | 50 digits |
| One-loop corrections (`virtual`, `loop_amplitude`) | FeynCalc + Package-X: e+e- -> mu+mu- in QED (vertices, vacuum polarization, IR-divergent boxes), H -> b b~ with the gluon loop, Z -> nu nu~ with all electroweak loops, the electron, gluon and ghost self-energies and the ghost-gluon vertex at xi = 1 and 3 (`tests/test_feyncalc_loop.sage`, `tests/test_diagrams.sage`) | poles and finite parts, 1e-15 |
| Textbook formulas | Peskin and Schroeder, read on the printed page: (5.11)-(5.13), (5.21), (5.87), (5.105)-(5.106), (7.90)-(7.91), (10.41), the QCD beta function (16.73)-(16.85) from feynsage's own counterterms, (21.107)-(21.108) (`tests/test_books.sage`, `tests/test_process.sage`) | exact |
| Diagram generation | FeynArts 3.12: topology counts (11 cases up to two loops) and diagram counts (13 processes up to 388 one-loop diagrams); FeynCalc 10.2 for e+e- -> mu+mu- in the SM and Bhabha; an explicit helicity sum for e+e- -> W+W- | identical counts; 20 digits |

Package-X and LoopTools both fail at one point we found (C0(-4,-3,-1; 0,0,1)). See
`tests/looptools/README.md`.

Run the tests with `sage tests/<name>.sage` from this folder.

More: the full description of every module, all the checks, speed and limits are in
[`docs/DETAILS.md`](docs/DETAILS.md).

## Licence

MIT, see `LICENSE`.
