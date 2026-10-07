# Dirac algebra with feynsage

This page shows how to do the Dirac, Lorentz and colour algebra of a QFT course with feynsage. You
write Python and Sage, in the notation of the textbook. FORM does the work in the background and the
answer comes back as an ordinary Sage expression. You never write FORM code.

```python
from sage.all import *
from feynsage import *

p, q = momenta("p q")
mu, nu = lorentz_indices("mu nu")

dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))
# 4*Comp(p, nu)*Comp(q, mu) + 4*Comp(p, mu)*Comp(q, nu) - 4*Dot(p, q)*Metric(mu, nu)
```

That is tr(p̸ γ^μ q̸ γ^ν) = 4(p^μ q^ν + p^ν q^μ - g^{μν} p·q). In a notebook `show(...)` prints it
with γ^μ, p̸, p·q and g^{μν}.

## Why FORM is used

FORM (J. Vermaseren) is a program made for very long expressions of exactly this kind. A trace of 12
slashed momenta has 10395 terms in d dimensions and FORM does it in about 14 ms. feynsage only builds
the problem from your Python objects, sends it to FORM and reads the answer back into Sage. The
algebra itself is never redone in Python.

One FORM process is started the first time and kept running, so a small trace costs about 0.1 ms.
Starting a new FORM process for each call would cost about 7.6 ms (measured on an M2 laptop, see the
table at the end).

## Momenta and Lorentz indices

```python
p, q, k = momenta("p q k")
mu, nu, rho, sigma_ = lorentz_indices("mu nu rho sigma_")
```

Both are plain Sage symbols, so `p - k`, `2*p` and `x*p + (1 - x)*q` work as usual. feynsage
remembers which symbols are momenta and which are indices. A name can be only one of the two.

`sigma_` has an underscore because `sigma` is also the function σ^{μν} below. Any name that Python
accepts is fine, FORM never sees your names directly.

## Gamma matrices and slashed momenta

| feynsage | textbook |
|---|---|
| `gamma(mu)` | γ^μ |
| `slash(p)`, `slash(p - k)` | p̸, (p̸ - k̸) |
| `gamma5()` | γ^5 |
| `PL()`, `PR()` | (1 - γ^5)/2, (1 + γ^5)/2 |
| `sigma(mu, nu)` | σ^{μν} = (i/2)[γ^μ, γ^ν] |
| `one()` | the unit matrix |

They multiply in the order you write them. `gamma(mu) * gamma(nu)` is never turned into
`gamma(nu) * gamma(mu)`. Numbers and Sage symbols mix freely, and a mass is added to a slash as in
the book.

```python
m = var('m')
slash(p) * gamma(mu) * (slash(q) + m)
# slash(p)*gamma(mu)*(slash(q) + m)
latex(_)
# {\not{p}} \gamma^{\mu} \left({\not{q}} + m\right)
```

Write `slash(p) + m`, not `slash(p + m)`. The second one raises an error that says so. `gamma(5)`
with a number is still Sage's Euler Γ, so `gamma(5)` is 24 and old Sage code keeps working.

## Dirac traces

```python
dirac_trace(gamma(mu) * gamma(nu))
# 4*Metric(mu, nu)
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sigma_))
# 4*Metric(mu, sigma_)*Metric(nu, rho) - 4*Metric(mu, rho)*Metric(nu, sigma_) + 4*Metric(mu, nu)*Metric(rho, sigma_)
dirac_trace((slash(p) + m) * gamma(mu) * (slash(q) + m) * gamma(nu))
# 4*m^2*Metric(mu, nu) + 4*Comp(p, nu)*Comp(q, mu) + 4*Comp(p, mu)*Comp(q, nu) - 4*Dot(p, q)*Metric(mu, nu)
```

`expr.trace()` is the same as `dirac_trace(expr)`. The trace of the unit matrix is 4.

Every Lorentz index that appears twice is summed, also inside the trace.
`dirac_trace(gamma(mu) * gamma(mu))` is 16 in four dimensions. A list is traced in one FORM run and
gives a list, which is the fastest way to do many traces:

```python
dirac_trace([gamma(mu) * gamma(mu), slash(p) * slash(p), 7])
# [16, 4*Dot(p, p), 28]
```

## Dimension

The default is four dimensions. `dim=d` (any symbol, or the string `"d"`) works in d dimensions.

```python
d = var('d')
dirac_trace(gamma(mu) * gamma(nu) * gamma(mu) * gamma(nu), dim=d)
# -4*d^2 + 8*d
```

## What comes back

The answer is a Sage expression made of these functions.

| printed | meaning | make it yourself with |
|---|---|---|
| `Dot(p, q)` | p·q (`Dot(p, p)` is p^2) | `dot(p, q)` |
| `Comp(p, mu)` | p^μ | `comp(p, mu)` |
| `Metric(mu, nu)` | g^{μν} | `metric(mu, nu)` |
| `Epsilon(mu, nu, p, q)` | ε^{μνρσ} with p and q contracted into two slots | `epsilon(mu, nu, p, q)` |

The printed names also work as input, so you can copy an output and paste it back. `dot` is
symmetric, `epsilon` is antisymmetric and both are linear in momenta, so
`dot(p + k, q)` is `Dot(k, q) + Dot(p, q)`.

Upper and lower indices are not told apart. A repeated index means g_{μν} is put in, as in the book.

## Contracting indices and kinematics

`contract` sums repeated indices of any such expression. It also reads the output of the older
string function `form.compute`.

```python
contract(metric(mu, nu) * comp(p, mu) * comp(q, nu))
# Dot(p, q)
contract(dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu)) * comp(k, mu) * comp(k, nu))
# 8*Dot(k, p)*Dot(k, q) - 4*Dot(k, k)*Dot(p, q)
contract(epsilon(mu, nu, rho, sigma_) * epsilon(mu, nu, rho, sigma_))
# -24
```

Scalar products are ordinary Sage objects, so kinematics is put in with `.subs`:

```python
s = var('s')
r = dirac_trace(slash(p) * slash(q) * slash(p) * slash(q))
# 8*Dot(p, q)^2 - 4*Dot(p, p)*Dot(q, q)
r.subs({dot(p, p): m^2, dot(q, q): 0, dot(p, q): s/2})
# 2*s^2
```

For a big expression it is faster to let FORM put the kinematics in before the answer comes back to
Sage. Give the same dictionary as `rules`. A rule can also replace a momentum, for momentum
conservation:

```python
dirac_trace(slash(p) * slash(q) * slash(p) * slash(q), rules={dot(p, p): m^2, dot(q, q): 0, dot(p, q): s/2})
# 2*s^2
p1, p2, k1, k2 = momenta("p1 p2 k1 k2")
dirac_trace(slash(k2) * slash(p1), rules={k2: p1 + p2 - k1})
# -4*Dot(k1, p1) + 4*Dot(p1, p1) + 4*Dot(p1, p2)
```

Momentum rules are used before the trace and scalar-product rules after it.

## Open lines

`contract` on a product of Dirac matrices puts the metric and the momenta into the gamma matrices and
then simplifies the line. Repeated indices are summed and p̸ p̸ = p^2.

```python
contract(gamma(mu) * gamma(nu) * gamma(mu))
# -2*gamma(nu)
contract(gamma(mu) * gamma(nu) * gamma(mu), dim=d)
# (-d + 2)*gamma(nu)
contract(metric(mu, nu) * gamma(mu) * gamma(nu))
# 4
```

FORM itself does not simplify open lines. This step uses the ordering of `feynsage.dirac`, which is
written in the conventions of Package-X.

## Gamma5

The conventions are those of Peskin and Schroeder: γ^5 = iγ^0γ^1γ^2γ^3 (their eq. 3.68),
tr(γ^μγ^νγ^ργ^σγ^5) = -4iε^{μνρσ} (eq. 5.5) and ε^{αβγδ}ε_{αβγδ} = -24 (eq. 5.6). Together they mean
ε^{0123} = +1. The tests check all 256 index values of the trace against explicit 4×4 matrices.

```python
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sigma_) * gamma5())
# -4*I*Epsilon(mu, nu, rho, sigma_)
dirac_trace(slash(p) * gamma(mu) * gamma5() * slash(q) * gamma(nu))
# 4*I*Epsilon(mu, nu, p, q)
```

In four dimensions FORM traces γ^5 exactly. In d dimensions there is no unique γ^5 and feynsage does
not choose one quietly. It uses an anticommuting γ^5 (NDR), the same as `feynsage.dirac` and
Package-X, only where the answer does not depend on the choice: an even number of γ^5 drops out and an
odd number with fewer than four other gamma matrices gives zero. Anything else is refused.

```python
dirac_trace(PL() * gamma(mu) * PR() * gamma(nu), dim=d)
# 2*Metric(mu, nu)
dirac_trace(gamma(mu) * gamma(nu) * gamma(rho) * gamma(sigma_) * gamma5(), dim=d)
# DiracError: a trace with gamma5 and 4 other gamma matrices in d dimensions is ambiguous with an
# anticommuting gamma5 (NDR).  Use dim=4, or do the gamma5 part in 4 dimensions.
```

The option `gamma5_scheme="NDR"` names this choice. Other schemes are not written yet and asking
for one raises an error. In the same way an ε tensor contracted in d dimensions needs `eps_in_d=True`,
because that is also a scheme choice.

## Spinors and fermion lines

`u(p, m)` and `v(p, m)` are columns and stand at the right end of a line. `ubar(p, m)` and
`vbar(p, m)` are rows and stand at the left end. A line with a spinor at both ends is a number.

```python
p1, p2 = momenta("p1 p2")
J = ubar(p2, m) * gamma(mu) * u(p1, m)
conjugate(J)
# ubar(p1, m)*gamma(mu)*u(p2, m)
```

`conjugate` reverses the matrices and takes γ^0 Γ^† γ^0 of each, so γ^5 changes sign, i becomes -i,
P_L and P_R are exchanged and σ^{μν} stays. Symbols are taken as real. List complex couplings in
`complex_symbols=`. An index that is summed inside the amplitude gets a new name in the conjugate
(`mu` becomes `mu_c`), so that M·M* never has one index four times.

`contract` uses the Dirac equation at the spinors, p̸u(p) = m u(p), ū(p)p̸ = m ū(p) and the same
with -m for v:

```python
contract(ubar(p2, m) * slash(p2) * gamma(mu) * slash(p1) * u(p1, m))
# m^2*ubar(p2, m)*gamma(mu)*u(p1, m)
```

Products that make no sense are stopped with an error in plain words, for example
`gamma(mu) * ubar(p, m)` gives "ubar(p, m) is a row spinor; it must stand at the left end of a line".

Two or more fermion lines are just multiplied. feynsage numbers the FORM lines itself. A relative
minus sign between diagrams (identical fermions) comes from Wick's theorem and you put it in M
yourself, feynsage never guesses it.

## Spin sums: e+ e- → μ+ μ-

`spin_sum` uses Σ u ū = p̸ + m and Σ v v̄ = p̸ - m, joins the bilinears of each term into closed
loops and lets FORM do one trace per loop.

```python
p1, p2, k1, k2 = momenta("p1 p2 k1 k2")
M_ = var('M_')
M = vbar(p2) * gamma(mu) * u(p1) * ubar(k1, M_) * gamma(mu) * v(k2, M_)
spin_sum(M * conjugate(M)) / 4
# 8*M_^2*Dot(p1, p2) + 8*Dot(k1, p2)*Dot(k2, p1) + 8*Dot(k1, p1)*Dot(k2, p2)
```

With p = p1, p' = p2, k = k1, k' = k2 this is Peskin and Schroeder eq. (5.10) without its factor
e^4/q^4. The average over initial spins (the 1/4 here) is left to you.

## Photons and polarization sums: Compton scattering

```python
kp, pp = momenta("kp pp")
e_in = polarization("e_in", k)        # ε(k)
e_out = polarization("e_out", kp)     # ε(k'); its conjugate ε*(k') is the vector e_out_c
eo = momenta("e_out_c")
Mc = (ubar(pp, m) * slash(eo) * (slash(p) + slash(k) + m) * slash(e_in) * u(p, m) / (2*dot(p, k))
      + ubar(pp, m) * slash(e_in) * (slash(p) - slash(kp) + m) * slash(eo) * u(p, m) / (-2*dot(p, kp)))
r = spin_sum(Mc * conjugate(Mc), rules={pp: p + k - kp})
r = polarization_sum(r, [e_in, e_out]) / 4
```

After p^2 = m^2, k^2 = k'^2 = 0 and k·k' = p·k - p·k' this is Peskin and Schroeder eq. (5.87) with
e = 1. The default is Σ ε^μ ε^{*ν} → -g^{μν}. `gauge="physical", n=p` uses the sum over the two
physical polarizations with a reference vector n and gives the same answer (the Ward identity). A
massive vector, `polarization("e", k, mass=MW)`, uses -g^{μν} + k^μ k^ν / M^2.

## Colour

```python
a, b, c = color_indices("a b c")       # gluon indices
i, j = quark_colors("i j")             # quark indices
color_factor(color_chain(a, a, i=i, j=j))
# 1/2*N*Kron(i, j) - 1/2*Kron(i, j)/N              (T^a T^a)_ij = C_F δ_ij
color_factor(f_color(a, b, c) * f_color(a, b, c), N=3)
# 24
```

`T_color(a, i, j)` is (T^a)_{ij} with Tr(T^a T^b) = δ^{ab}/2 and `f_color(a, b, c)` is f^{abc}.
`color_trace(a, b, c)` is Tr(T^a T^b T^c). FORM reduces everything with the SU(N) Fierz identity, so the
answer is in N. Colour factors may sit inside amplitudes. `spin_sum` and `contract` sum them along
with the Lorentz indices and `conjugate` turns (T^a)_{ij} into (T^a)_{ji}. The tests compare with
explicit Gell-Mann matrices.

## Conventions: standard or the note, Minkowski or Euclidean

```python
set_convention("standard")       # Peskin and Schroeder, the dimension symbol is d
set_convention("note")           # the lecture note, D = 4 - 2 eps
spacetime_dimension()            # d or D, by the convention
set_convention(euclidean=True)   # Euclidean algebra
```

Both conventions use the metric (+, -, -, -). With `euclidean=True` the algebra is
{γ_μ, γ_ν} = 2δ_{μν}, γ^5 = γ_1γ_2γ_3γ_4, tr(γ_μγ_νγ_ργ_σγ^5) = 4ε_{μνρσ} with ε_{1234} = +1 and
ε·ε = +24, checked with explicit Euclidean matrices. The metric prints as `Delta(mu, nu)` and LaTeX
writes squares as p_E^2, as in the note. Mixing a Minkowski metric into a Euclidean calculation is an
error. Spin sums are written for Minkowski spinors only.

To move a scalar result between the two, use the rule of the note, p·q = -p_E·q_E (so s = -p_E^2):

```python
to_euclidean(dot(p, p) - m^2)
# -Dot(p, p) - m^2
```

## When something goes wrong

Errors talk about physics, not FORM. Some you may meet:
"gamma(muu) is Euler's Gamma function of the symbol muu here, because muu is not a Lorentz index",
"the index mu appears more than twice in one term" and
"spin_sum: u(p1) has no partner ubar(p1) to be summed with".

`debug=True` prints the FORM program and FORM's own output:

```python
dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(mu), debug=True)
# --- FORM program ---
# Vectors p,q;
# Indices mu;
# Off statistics;
# Local FSE1 = g_(1,p,mu,q,mu);
# trace4,1;
# contract;
# Format nospaces;
# Print;
# .end
# --- FORM output ---
# FSE1=-8*p.q;
```

If FORM itself fails you get a `FormError` with the program and FORM's message attached
(`err.program`, `err.output`).

## Speed

`benchmarks/qft_vs_form.sage` traces n slashed momenta and compares with FORM run by hand on the
same program. On an M2 laptop:

| n | dim | terms | FORM by hand (s) | new FORM process per call (s) | kept FORM process (s) |
|---|---|---|---|---|---|
| 8 | 4 | 105 | 0.007 | 0.008 | 0.001 |
| 10 | d | 945 | 0.008 | 0.015 | 0.007 |
| 12 | 4 | 4383 | 0.011 | 0.041 | 0.034 |
| 12 | d | 10395 | 0.014 | 0.095 | 0.082 |

Up to about a thousand terms feynsage takes no more time than starting FORM by hand, because the FORM
process is kept. Above that most of the time goes into building the Sage expression (4383 terms take
0.034 s against 0.011 s for FORM alone). If you only need
a number at the end, put the kinematics in with `rules=` so that FORM returns a short answer.
`set_backend("subprocess")` (or the environment variable `FEYNSAGE_BACKEND=subprocess`) starts a new
FORM process for every call instead. `threads=8` runs FORM's multi-threaded tform on a big job.

## Writing FORM yourself

The old interface stays. `form.compute(...)`, `form.dirac_trace(...)` and `form.run_form(code)` take
strings, and `run_form` runs any FORM program you write. Their results can be passed to `contract`.

## How it is checked

`sage tests/test_qft_dirac.sage` runs 90 checks in about 4 seconds. It compares traces with
explicit Dirac matrices in Minkowski and Euclidean space, with FORM programs run directly through
`form.compute`, with Peskin and Schroeder eqs. (5.5), (5.6), (5.9), (5.10) and (5.87), and with
Gell-Mann matrices for colour. It also checks that the wrong inputs listed above are refused.
