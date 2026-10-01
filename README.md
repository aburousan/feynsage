# feynsage

Feynman integrals in SageMath, written from scratch.

- **One loop, like Package-X.** Tensor integrals reduced to A0, B0, C0, D0, exactly, in the same
  conventions. C0 and D0 in closed form (logarithms and dilogarithms). Soft and collinear poles
  written out. Works also when the Gram determinant is zero.
- **Many loops, like LiteRed, FIRE and Kira.** Integral families from strings, IBP reduction
  (exact, or with finite fields), symmetries found by itself, master integrals.
- **Graphs.** The Symanzik polynomials U and F from spanning trees and 2-forests, checked against
  two other methods. Pictures of diagrams, trees and sectors.
- **Dirac algebra with FORM**, gamma_5 included.
- Every function can explain its own output: `explain=True` or `info(f)`.

I wrote feynsage in 2024 for my MSc project. After the talk "Feynman integrals at one loop" by
Prof. B. Ananthanarayan (NISER, 2026) it was refined, checked against the lecture notes of that
talk and published.
A short talk with real computations is in [`slides/feynsage_talk.pdf`](slides/feynsage_talk.pdf).

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

## Notebooks

- `tutorials/`: five short step-by-step notebooks for a first contact, each with exercises and
  solutions: first steps (A0, B0, numbers, plots), one-loop tensors, graph polynomials at one, two
  and three loops, families, IBP and master integrals (the kite as 6 zeta(3)), Dirac traces with FORM.
  Start with [`tutorials/README.md`](tutorials/README.md).
- `examples/feynsage_walkthrough.ipynb`: the whole package step by step, with a picture for every
  step: graph polynomials (up to a three-loop box), IBP reduction (finite fields and fast exact),
  one-loop integrals with IR poles, a real prediction (the neutral pion lifetime, 8.35e-17 s against
  the measured 8.43e-17 s) and the QED of Chluba's thesis (Compton, double Compton, the infrared
  divergence and its cancellation), with Feynman diagrams.
- `examples/peskin_examples.ipynb`: vacuum polarisation and the electron g-2 from Peskin and
  Schroeder, exactly.
- `examples/chluba/`: Compton and double Compton scattering from J. Chluba's thesis (CMB spectral
  distortions), with feynsage and with FeynCalc, and the cancellation of the infrared divergence in
  dimensional regularisation.

## Is it right?

The tested results are checked against another program, a known answer or a direct integral:

| Part | Checked against | Agreement |
|---|---|---|
| A0, B0, C and D tensors | Package-X 2.1.1, 190 points | every digit of its values |
| IR-divergent C0, D0 | Package-X, 118 points; a photon mass where Package-X has no number | double precision |
| D0, also with massless lines | Package-X (146 points), LoopTools, direct integration | 1e-9 (Package-X's own limit), 1e-15 against the others |
| IBP reduction | Kira 3.1, LiteRed | identical |
| U and F | sir's notebook (uf-new-short.nb); FeynCalc | identical |
| QED | Peskin and Schroeder | Ward identity and F2(0) = alpha/2pi exactly |

Package-X and LoopTools both fail at one point we found (C0(-4,-3,-1; 0,0,1)). See
`tests/looptools/README.md`.

Run the tests with `sage tests/<name>.sage` from this folder.

More: the full description of every module, all the checks, speed and limits are in
[`docs/DETAILS.md`](docs/DETAILS.md).

## Licence

MIT, see `LICENSE`.
