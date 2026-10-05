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
- **Dirac algebra with FORM**, gamma_5 included.
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

- `examples/feynsage_first_steps.ipynb`: for a first course in QFT. The recipe of every Feynman calculation and
  one short, explained step for each part: a tree diagram ($e^+e^- \to \mu^+\mu^-$), the tadpole and the bubble,
  the running coupling of $\lambda\varphi^4$, tensor integrals, triangles and boxes, the electron $g - 2$,
  Feynman parameters and a two-loop reduction. Every result is checked against Peskin and Schroeder or a second method.
- [`docs/tutorial/feynsage_tutorial.pdf`](docs/tutorial/feynsage_tutorial.pdf): the step-by-step tutorial, in the
  style of the lecture note. It follows the topics of the talk (the tadpole, IBP, the bubble and its
  differential equation solved with its boundary condition, Feynman parameters, graph polynomials) and
  ends with a real prediction, the decay rate of the neutral pion into two photons. Every step has its code, its real output and
  an explanation of each function used. Reductions can be drawn as equations of diagrams (`r.draw`).
- `examples/feynsage_walkthrough.ipynb`: the whole package step by step, with a picture for every
  step: graph polynomials (up to a three-loop box), IBP reduction (finite fields and fast exact)
  with reductions drawn as diagrams, Dirac traces from Sage notation with FORM, the bubble
  differential equation solved from its boundary condition, one-loop integrals with IR poles, a real
  prediction (the decay rate of the neutral pion into two photons, 7.79 eV) and the QED of
  Chluba's thesis (Compton, double Compton, the infrared divergence and its cancellation), with
  Feynman diagrams.
- `examples/peskin_examples.ipynb`: vacuum polarisation and the electron g-2 from Peskin and
  Schroeder, exactly.
- `examples/chluba/`: Compton and double Compton scattering from J. Chluba's thesis (CMB spectral
  distortions) with feynsage and with FeynCalc. Also the cancellation of the infrared divergence in
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

Package-X and LoopTools both fail at one point we found (C0(-4,-3,-1; 0,0,1)). See
`tests/looptools/README.md`.

Run the tests with `sage tests/<name>.sage` from this folder.

More: the full description of every module, all the checks, speed and limits are in
[`docs/DETAILS.md`](docs/DETAILS.md).

## Licence

MIT, see `LICENSE`.
