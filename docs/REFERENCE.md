# feynsage function reference

Every public function and class of feynsage, module by module, with its signature and what it returns.
`from feynsage import *` loads the names marked **(top)**. The others are reached through their
module, for example `feynsage.oneloop.tadpole` or `feynsage.form.compute`.

Every function explains itself in more detail: `info(f)` prints what it computes and what its output
means and the main ones also take `explain=True`. Worked examples are in the tutorial
([`tutorial/feynsage_tutorial.pdf`](tutorial/feynsage_tutorial.pdf)) and in
[`../examples/feynsage_walkthrough.ipynb`](../examples/feynsage_walkthrough.ipynb). What each module
does and how it is checked is in [`DETAILS.md`](DETAILS.md).

Two normalisations are used.
- `oneloop`, `family`, `graph` and `ibp_reduce`: Euclidean, measure d^D l / pi^(D/2), propagators l^2 + m^2
  (`euclidean=False` in `family` and `graph` gives Minkowski propagators instead).
- `pv`, `scalar`, `ir`: Package-X and LoopTools, mu^(2 eps) e^(eps gamma_E) Int d^D l / (i pi^(D/2)) in Minkowski space.

Contents: [easy](#easy-families-graphs-and-reduction) ·
[oneloop](#oneloop-closed-formulas) ·
[pv](#pv-one-loop-integrals-in-the-package-x-normalisation) ·
[dirac](#dirac-fermion-lines-projectors-package-xs-fermionline-layer) ·
[scalar](#scalar-c0-and-d0-in-closed-form) ·
[ir](#ir-infrared-divergent-triangles-and-boxes) ·
[form](#form-dirac-algebra-with-form) ·
[plotting](#plotting-figures-and-diagrams) ·
[graph](#graph-feynmangraph) ·
[family](#family-integralfamily) ·
[laporta](#laporta-reducer) ·
[ff](#ff-finite-field-reduction) ·
[tex, explain, momenta](#tex-explain-and-momenta)

---

## easy: families, graphs and reduction

The short front ends. Everything is written as plain strings.

| Function | What it does |
|---|---|
| `family(props, kin=None, loops=None, euclidean=False, name="F", explain=False)` **(top)** | An integral family from a list of propagators. A propagator is a momentum `"l + p"` or a pair `("l + p", "m")` for a massive line. `kin` gives the external scalar products, for example `{"p^2": "s"}`. A Sage name cannot contain `^`, so write a value like k^2 as a symbol such as `"kk"`. `loops` names the loop momenta (by default every name that is not external). Returns a `Family`. |
| `graph(edges, legs, kin=None, euclidean=False, explain=False)` **(top)** | A Feynman graph. `edges` is `"A-B, A-B:m"` (`:m` gives the line a mass), `legs` the incoming momentum at each vertex, for example `{"A": "p", "B": "-p"}`. Lines are numbered x1, x2, ... in the order written. Returns a `FeynmanGraph`. |
| `diagram(name=None)` **(top)** | A ready-made `FeynmanGraph`: `tadpole`, `bubble`, `bubble_mass`, `triangle`, `box`, `sunset`, `kite`, `phi4_two`, `vertex2`, `banana3`, `ladder3`, `triplebox`. `diagram()` prints the list. |
| `ibp_reduce(fam, targets, method="auto", symmetries_="auto", nproc=1, explain=False, verbose=False)` **(top)** | Reduces the target integrals to master integrals. `targets` is a list like `["F(2,1)", (1,2)]`, the integral named by its powers (the letter is the family name). `method` is `"ff"` (finite fields, fast, at most one symbolic invariant besides d), `"trimmed"` (exact, only the equations the targets need) or `"exact"` (exact on every seeded equation). `"auto"` chooses. `nproc > 1` samples in parallel. Returns a `Reduction`. |
| `symmetries(fam)` **(top)** | The permutations of the propagators that leave U + F unchanged (Pak's criterion), as generators. `ibp_reduce` calls it by itself. |
| `kinematics(externals, kin=None, euclidean=False, extra=())` **(top)** | A `Kinematics` object from `{"p^2": "s", "p1.p2": "t/2"}`. A product that is not given becomes a new symbol. |

**`Family`**: an `IntegralFamily` (see [family](#family-integralfamily)) that also remembers how it was written.
`.info()` prints its propagators and conventions.

**`Reduction`**: what `ibp_reduce` returns.

| Member | What it is |
|---|---|
| `r[target]` | `{master: coefficient}` for one target, given as `"J(2,1)"` or `(2, 1)` |
| `r.table` | `{target: {master: coefficient}}` for all targets |
| `r.masters` | the list of master integrals |
| `r.method`, `r.seconds` | the method used and the time it took |
| `r.info()` | a summary |
| `r.draw(graph, targets=None, rename=None, size=1.5)` | the reduction as equations of diagrams. `targets` as `"J(2,1)"` or `(2, 1)`. `graph` must have its lines in the order of the propagators. A line with power 0 is shrunk to a point, each extra power is a dot, `rename={"kk": "k^2"}` prints a variable under another name. Returns a matplotlib figure. |

In a notebook a `Reduction` prints as typeset equations.

---

## oneloop: closed formulas

Euclidean, measure d^D l / pi^(D/2). `D` is the symbol of the dimension (`from feynsage.oneloop import D`).

| Function | What it returns |
|---|---|
| `tadpole(n, m2)` | Int 1/(l^2 + m^2)^n = Gamma(n - D/2) (m^2)^(D/2 - n) / Gamma(n). `m2` is the mass squared. |
| `bubble_massless(a, b, p2)` | the massless bubble with powers a and b at p^2 = `p2`, G(a, b) (p^2)^(D/2 - a - b) |
| `G(a, b)` | the G-function of the massless bubble, a ratio of Gamma functions |
| `bubble_equal_mass(p2, m2)` | the equal-mass bubble, Gamma(2 - D/2) m^(D-4) 2F1(2 - D/2, 1; 3/2; -p2/(4 m2)) |
| `bubble_one_mass(a, b, q2, m2)` | the bubble with one massless line (power a) and one massive line (power b), as a 2F1 |
| `triangle_onshell(a1, a2, a3, Q2)` | the massless triangle with two on-shell legs and powers a1, a2, a3 |
| `expand_eps(expr, order=0, loops=1)` | the Laurent expansion at D = 4 - 2 eps up to eps^order, with the factor e^(L eps gamma_E). Gamma functions with a pole are made regular first. |

---

## pv: one-loop integrals in the Package-X normalisation

| Function | What it returns |
|---|---|
| `A0(m)` **(top)** | the tadpole, m^2 (1/eps + log(mu^2/m^2) + 1). `A0(0) = 0`. |
| `B0(s, m1, m2)` **(top)** | the two-point function at p^2 = s: 1/eps plus the finite part in closed form, with `LogM` and `DiscB` |
| `C0(s1, s12, s2, m0, m1, m2)` **(top)** | the scalar triangle. Propagators (l, m0), (l + p1, m1), (l + p2, m2), s1 = p1^2, s2 = p2^2, s12 = (p1 - p2)^2. IR finite: the symbol `C0(...)`, whose `.n()` gives more than 30 digits. IR divergent: the poles and the finite part. |
| `D0(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3)` **(top)** | the scalar box in the Package-X / LoopTools order, the same behaviour as `C0` |
| `PVA(r, m)` **(top)** | the tadpole coefficient as Package-X's `PVA[r, m]`: `PVA(0, m)` = A0, `PVA(1, m)` = A00 |
| `PVB(r, n, s, m1, m2, explain=False, series=None, disc=None)` **(top)** | the Passarino-Veltman coefficient as Package-X's `PVB[r, n, ...]`: `PVB(0,0)` = B0, `PVB(0,1)` = B1, `PVB(1,0)` = B00. `series=(x, x0, n)` gives the Taylor series to order (x - x0)^n (Package-X's LoopRefineSeries), `disc=s` the discontinuity across the cut in s (Package-X's Part -> Discontinuity[s]). The same two options exist for `PVC` and `PVD`; `PVD` also takes `disc=(s, t)` for the double spectral function. |
| `PVC(r, n1, n2, s1, s12, s2, m0, m1, m2, explain=False)` **(top)** | the same for triangles, as `PVC[r, n1, n2, ...]` |
| `PVD(r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, explain=False)` **(top)** | the same for boxes, as `PVD[...]` |
| `loop(numerator, *props, kin=None, euclidean=False, explain=False)` **(top)** | a one-loop tensor integral reduced to A0, B0, C0, D0. A propagator is `["l + p", "m"]`; writing the same propagator twice raises its power (reduced by integration by parts). Up to five propagators: a pentagon is written through boxes (exact up to O(eps), tensors up to rank 5). The numerator may contain `l^2`, `l.p`, `l^mu`, `g(mu,nu)` and products. Returns a `LoopResult`. |
| `loop_diff(expr, x)` **(top)** | d expr/dx for an expression with A0, B0, C0, D0, written again in the same functions; the derivatives of C0 and D0 are exact (integration by parts) |
| `loop_series(expr, (x, x0, n), ...)` **(top)** | the Taylor series of a result about a point where its closed form is regular. About zero external momenta use the `series=` option of `PVB`, `PVC`, `PVD`, which expands the Feynman-parameter integral. |
| `explicit(expr)` **(top)** | every `C0(...)` and `D0(...)` with numerical arguments replaced by its closed form in dilogarithms |
| `finite_part(expr, mu_value=None)` **(top)** | the eps^0 part. `mu_value` sets the scale mu to a number. |
| `uv_part(expr)` **(top)** | the coefficient of 1/eps |
| `pole_parts(expr)` **(top)** | (coefficient of 1/eps^2, coefficient of 1/eps) |
| `c0_numeric(s1, s12, s2, m0, m1, m2, n=260, tol=1e-06)` **(top)** | C0 by numerical integration on two grids, as a check |
| `d0_numeric(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, n=90, tol=1e-05)` **(top)** | D0 the same way |
| `LogM(z)`, `DiscB(s, m1, m2)` **(top)** | the functions in the finite part of `B0`, as in Package-X: log(z - i0) and the discontinuity function of the bubble |
| `eps`, `mu` **(top)** | the symbols eps (D = 4 - 2 eps) and mu (the scale) |

**`LoopResult`**: what `loop` returns. `.coefficients()` gives the coefficient of each tensor structure,
`.masters()` the same before the eps expansion (exact in d), `.discontinuity(s)` (or `(s, t)`) the cut, `.info()` a summary.

Kinematic helpers and branch-aware functions, with Package-X's meaning:

| Function | What it returns |
|---|---|
| `kallen(a, b, c)` **(top)** | the Kallen function a^2 + b^2 + c^2 - 2ab - 2ac - 2bc |
| `kibble(s1, s2, s3, s4, s12, s23)` **(top)** | the Kibble polynomial of a four-point function |
| `mandelstam(momenta, masses, stu)` **(top)** | the scalar products of p1 + p2 -> p3 + p4 in s, t, u, as a `kin` dictionary |
| `disc_expand(expr)` **(top)** | DiscB replaced by its definition with sqrt and log |
| `Ln(x, a)`, `DiLog(x, a)` **(top)** | log(x + i a 0) and Li2(x + i a 0) for real x, the side of the cut given by the sign of a |
| `continued_dilog(x1, a1, x2, a2)` **(top)** | the Beenakker-Denner continued dilogarithm, as a number |

`LogM`, `DiscB`, `Ln`, `DiLog` evaluate to any precision asked for with `.n(digits=...)`.

---

## dirac: fermion lines, projectors (Package-X's FermionLine layer)

`from feynsage import dirac`. Dirac matrices are lists of factors as in `form.dirac_trace`: an index `"mu"`,
slashed momenta plus a mass `"p - k + m"`, `"g5"`, `"PL"`, `"PR"`. d-dimensional algebra with an
anticommuting gamma_5. A spinor is `("u" or "v", momentum, mass)`. Results are dictionaries
`{(structure, chirality): coefficient}` with the structures `()` (the unit matrix), `(("i", mu),)` (gamma^mu),
`("sigma", a, b)` (sigma^{ab}, a and b an index or a momentum) and ordered products; the chirality is
`"1"`, `"5"`, `"L"` or `"R"`. In the coefficients `p__mu` is p^mu and `g__mu__nu` is g^{mu nu}.

| Function | What it does |
|---|---|
| `line_expand(left, factors, right, sp=None, gordon=True)` **(top)** | Package-X's FermionLineExpand: contractions, the Dirac equation at both ends, on-shell momenta, sigma basis and (by default) the Gordon identities |
| `loop_line(left, factors, right, *props, kin=None, gordon=True)` **(top)** | a one-loop integral over a fermion line (Package-X's LoopIntegrate + LoopRefine of a FermionLine). The factors may contain the loop momentum `l`; the line is simplified exactly in d before the eps expansion. Returns a `LineResult`. |
| `loop_matrix(factors, *props, kin=None)` **(top)** | the same without spinors (a self-energy matrix) |
| `projector(name, mu, p1, m1, p2, m2, q2="q2")` **(top)** | Package-X's Projector: vertex form factors F1 (Dirac), F2 (Pauli), F3, G1 (Anapole), G2 (EDM), G3, SachsElectric, SachsMagnetic, the chiral AL ... CR; with `mu=None` the densities S, P; with `p=..., m=...` the self-energy A, B, C, E and AL, BL, AR, BR. Built by linear algebra, exact in d. |
| `form_factor(factors, P, *props, kin=None)` **(top)** | Int Tr[M(l) P]/prod D_i: a form factor of a loop, the projection done before the integral |
| `spur(M, P, momenta, sp)`, `trace(factors, momenta, sp)` | traces of products and of sums of products |
| `transverse(T, v, vsq)`, `longitudinal(T, v, vsq)` **(top)** | Package-X's Transverse and Longitudinal for a rank-2 `loop()` result |
| `chisholm(result)` **(top)** | the four-dimensional Chisholm identity on products of three gamma matrices, + i eps^{a b c s} g_s g5 with eps^{0123} = +1 (Peskin and Schroeder, Package-X) |
| `to_chiral(result)`, `to_g5(result)` **(top)** | switch between the 1, g5 and PL, PR forms |

| `line_product(*lines, chisholm=True, gordon=True, on_shell=True)` **(top)** | Package-X's FermionLineProduct with FermionLineExpand: each line is `(left, factors, right)` as in `line_expand`; an index in two lines is summed. Lines are reduced in the d-dimensional antisymmetric basis and contracted in d dimensions; `chisholm=True` maps rank 3 and 4 to their four-dimensional duals (eps^{0123} = +1, eps.eps contracted with the d-dimensional metric). Exact at d = 4; away from d = 4 Package-X keeps different evanescent terms for some index orders. Returns a `LineProduct` {((structure, chiral), ...): coefficient}, summed indices renamed fs1, fs2, ... |

---

## scalar: C0 and D0 in closed form

| Function | What it returns |
|---|---|
| `c0_closed(s1, s12, s2, m0, m1, m2, check=True)` | C0 as an exact expression in `polylog(2, ...)` and `log` of algebraic numbers, every branch fixed by the +i0. `check=True` compares with the numerical value. |
| `c0_value(s1, s12, s2, m0, m1, m2, full=False)` | the high-precision complex value from the same construction |
| `d0_closed(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, check=True)` | D0 in closed form (16 dilogarithms with all masses nonzero, fewer when masses vanish). IR-finite boxes only. |
| `d0_value(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, full=False)` | the high-precision complex value (about 30 digits) |
| `c0_values(arglist, nproc=None)`, `d0_values(arglist, nproc=None)` **(top)** | many triangles or boxes at once, on every core |

---

## closed: C0 and D0 with symbols (Package-X's C0Expand and D0Expand)

| Function | What it returns |
|---|---|
| `c0_expand(s1, s12, s2, m0, m1, m2)` **(top)** | C0 as an explicit formula in the arguments (Denner's twelve dilogarithms). The i0 is written in as the symbol fs_eps; named intermediates (alpha, y0_i, x_i+-) come in `.defs`. A `Conditional`, valid where lambda(s1, s12, s2) > 0. |
| `d0_expand(s1, ..., s23, m0, ..., m3, pair=(1, 3))` **(top)** | D0 as an explicit formula: Denner's sixteen dilogarithms with all masses nonzero (needs a real r for the line pair `pair`; the default is Package-X's condition on s23), Denner and Dittmaier's formulas when masses vanish (all real kinematics) |
| `expand_c0d0(expr)` **(top)** | every C0(...) and D0(...) inside expr replaced (C0Expand and D0Expand together) |
| `Conditional` | `.value`, `.defs`, `.conditions`; `.n(point)` the value (fs_eps, fs_del -> 0 at 600 bits), `.n_many(points)` on every core, `.holds(point)`, `.inline()` one expression, `.subs(...)` |
| `SqrtL`, `Eta`, `EtaT`, `ReSign` | the small functions of the formulas: the limit square root, eta, Denner's eta-tilde, the sign of a real part |

---

## ir: infrared-divergent triangles and boxes

The six triangles and sixteen boxes of Ellis and Zanderighi (2008), in the Package-X normalisation.
`C0` and `D0` call this module by themselves when the integral is IR divergent.

| Function | What it returns |
|---|---|
| `classify_c0(s1, s12, s2, m0, m1, m2)` | `None` for an IR-finite C0, else the Ellis-Zanderighi triangle number and its relabelled arguments |
| `classify_d0(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3)` | the same for D0 |
| `ir_coefficients(kind, args, exact=True)` | (c_-2, c_-1, c_0) at mu = 1 for `kind` `"C"` or `"D"` |
| `ir_value(kind, args)` | the same three coefficients as complex numbers |
| `with_mu(c)` | c_-2/eps^2 + c_-1/eps + c_0 with the scale mu put back |

---

## form: Dirac algebra with FORM

FORM must be installed (`install.sh` does it). Sage-side notation: `g(mu,nu)` the metric, `dot(p,q)`
the scalar product, `comp(p,mu)` the component p^mu, `eps(a,b,c,d)` the Levi-Civita symbol.

| Function | What it returns |
|---|---|
| `compute(expr="1", vectors=(), lines=(), rules=None, dim=4, show_code=False, levi_civita="form")` | traces and index contractions without writing FORM code. `expr` is a product such as `"eps(mu,nu,rho,sigma)*k1(rho)*k2(sigma)"`, `lines` the fermion lines to trace (each a list of factors as in `dirac_trace`), `rules` the scalar products to put in, for example `{"k1.k1": 0}`. Repeated indices are summed. `show_code=True` prints the FORM program. `levi_civita="usual"` uses the epsilon of Peskin and Schroeder (eps^0123 = +1, Tr[g^mu g^nu g^rho g^sigma g5] = -4i eps^{mu nu rho sigma}), written `Eps`. The default `"form"` keeps FORM's `e_` = -i epsilon. Returns a Sage expression. |
| `dirac_trace(factors, vectors, dim=4, levi_civita="form")` | the trace of one product. A factor is `"g5"`, an index such as `"mu"`, or a linear expression such as `"q - x*k1 + m"` (slashed momenta plus a mass). `dim=4` allows gamma_5, a symbol such as `"D"` works in D dimensions. |
| `trace(indices, dim="D")` | the trace of gamma_{i1} ... gamma_{in} in `dim` dimensions |
| `to_loop(expr)` | a `dirac_trace` result as a numerator string for `loop` |
| `run_form(code)` | runs a FORM program written by hand and returns its output |

The blog post [FORM for pedestrian QFT](https://rousan.netlify.app/pages/physics/blogs/form_pedestrian_qft_blog/)
explains FORM itself.

---

## plotting: figures and diagrams

`from feynsage import plotting`. Every function returns a matplotlib figure, saved with `fig.savefig("name.svg")`.

| Function | What it draws |
|---|---|
| `quick_plot(exprs, var_range, labels=None, points=160, parts="auto", title=None, ylim=None, nproc=1, figsize=None, mu_value=1, explain=False)` | one or more expressions over `(variable, a, b)`, real part and (if not zero) imaginary part. Expressions with `DiscB`, `LogM`, `C0`, `D0` and 1/eps poles work directly. |
| `draw_graph(g, labels=True, dots=None, momenta=True, figsize=(2.8, 2.8), R=1.0, removed=(), ax=None, title=None, names=None)` | one diagram. Massive lines thick and blue, line i labelled x_i, `dots={i: n}` puts n dots on line i, `removed` dashes lines. `g.plot()` calls it. |
| `draw_panels(g, removed_sets, titles=None, ncols=5, size=1.9, **kw)` | copies of a diagram side by side, each with its own removed lines (all spanning trees or all 2-forests) |
| `draw_sectors(g, sectors, titles=None, ncols=5, size=1.9, **kw)` | sectors of a diagram, with the lines of index 0 shrunk to points |
| `draw_reduction(table, g, targets=None, rename=None, size=1.5, fontsize=12)` | a reduction table as equations of diagrams (`Reduction.draw` calls it) |
| `curves(funcs, xrange, labels=None, points=400, xlabel="$x$", ylabel=None, title=None, ylim=None, figsize=(4.2, 2.8))` | several real functions on one axis |
| `mb_plane(left, right, c, left_label=None, right_label=None, close=None, figsize=(5.0, 3.0))` | the complex plane of a Mellin-Barnes integral: left and right poles and the contour Re z = c |
| `set_theme()` | the style of the notes for every later matplotlib and Sage plot |

---

## graph: FeynmanGraph

Made by `graph(...)` or `diagram(...)`.

| Method | What it returns |
|---|---|
| `g.U()` | the first Symanzik polynomial, from spanning trees |
| `g.F()` | the second Symanzik polynomial, F_0 plus U times the masses |
| `g.F0()` | the massless part of F, from 2-forests |
| `g.U_kirchhoff()` | U from the matrix-tree theorem, an independent check |
| `g.spanning_trees()` | the removed-line sets of all spanning trees |
| `g.two_forests()` | (removed lines, one of the two components) for all spanning 2-forests |
| `g.contract(lines)` | the graph with these lines shrunk to points and the original numbers of the remaining lines |
| `g.family(name="fam", loop_names=None)` | momentum routing: (propagators, loop names), ready for `IntegralFamily` |
| `g.plot(**kw)` | the diagram (see `draw_graph`) |

---

## family: IntegralFamily

Made by `family(...)`.

| Method | What it returns |
|---|---|
| `fam.UF(sector=None, names="x")` | (U, F) by completing the square: U = det M, F = U (J - Q.M^-1.Q) in Euclidean signature. `sector` is a 0/1 tuple, by default all lines. |
| `fam.ibp(a)` | all L(L + E) IBP identities for the power vector `a`, each a dict `{powers: coefficient}` |
| `fam.is_zero_sector(sector)` | `True` if the sector integrals vanish (Lee's criterion) |
| `fam.dot(u, v)` | the scalar product of two momenta in terms of the loop scalar products and the kinematics |
| `fam.sp_in_dens(j)` | scalar product number j written through the propagators |

---

## laporta: Reducer

The Laporta reducer behind `ibp_reduce`, for users who want to control the seeds.

`Reducer(family, symmetries=(), sector_symmetries=True, numerator_relations="auto", top=None)`

| Method | What it does |
|---|---|
| `run(rmax, smax=0, verbose=False)` | writes and solves the identities for all seeds with at most `rmax` dots and `smax` numerator powers |
| `reduce(a)` | the reduction of the integral with powers `a` |
| `masters(targets=None)` | the master integrals (of the targets, or of the whole seeded range) |
| `seeds(rmax, smax)` | the seed integrals, sector by sector |
| `set_top(integrals)` | makes the union of the sectors of these integrals the top sector |
| `canon(a)` | the representative of `a` under the symmetries |
| `sector(a)`, `is_zero(a)` | the sector of `a` and whether it vanishes |
| `symmetry_relations(a, polynomial=False)` | extra equations for an integral with numerators from momentum shifts |
| `relations_needed(masters)` | whether those extra equations are needed |

`weight(a)` is the ordering of integrals used to pick masters.

---

## ff: finite-field reduction

| Function | What it does |
|---|---|
| `reduce_ff(reducer, targets, rmax, smax=0, point=None, verbose=False, nproc=1)` **(top)** | the reduction done modulo large primes and rebuilt by Thiele interpolation and rational reconstruction, the way FIRE and Kira do it. Returns `{target: {master: coefficient}}`. |
| `reduce_exact_trimmed(reducer, targets, rmax, smax=0, verbose=False)` | exact Laporta on only the equations the targets need (found by one modular probe) |
| `ibp_templates(fam)` | the IBP identities with symbolic powers, derived once |
| `IBPSystem(reducer, rmax, smax=0)` | all identities of the seeds with polynomial coefficients. `.sample(targets, p, point)` solves them modulo p at one point. |

---

## parallel: every core

| Function | What it returns |
|---|---|
| `pmap(f, items, nproc=None)` **(top)** | `[f(x) for x in items]` on worker processes (fork); f may be a lambda or a closure. Runs serially when the items are too few or too fast to pay for the workers. |
| `cores()`, `set_nproc(n)` **(top)** | the number of workers (set_nproc, else FEYNSAGE_NPROC, else all cores); `set_nproc(1)` switches parallel work off |
| `scan(expr, var, values, nproc=None)` **(top)** | the values of an expression (with DiscB, LogM, C0, D0) at many points |
| `loop_many(jobs, nproc=None)` **(top)** | many `loop()` calls: jobs are (numerator, propagators, options) |

Used by default (nproc=None): `ibp_reduce`/`reduce_ff` sample points, `quick_plot`, `scan`, `c0_values`, `d0_values`, `Conditional.n_many`, `loop_many`. FORM: `form.run_form(code, threads=n)` runs `tform -w n` (default from FEYNSAGE_FORM_THREADS, else 1).

---

## tex, explain and momenta

| Function | What it does |
|---|---|
| `eq(lhs, expr, per_line=2)` **(top)** | shows `lhs = expr` as a typeset equation in Jupyter |
| `tex.laurent(expr, per_line=2)` | LaTeX of c_-2/eps^2 + c_-1/eps + c_0, poles first, finite part collected by functions |
| `tex.lines(expr, per_line=3)` | LaTeX of a long sum broken into rows |
| `tex.grouped(expr)` | the terms of `expr` collected by their transcendental factor |
| `info(obj)` **(top)** | prints the description of a feynsage function, class or result |
| `explain.explainable(f)` | gives a function the `explain=True` keyword |
| `mom(**kw)` **(top)** | a momentum as a linear combination: `mom(l1=1, p1=-1)` is l1 - p1 |
| `Kinematics(externals, invariants, rules, euclidean=False, extra_vars=())` **(top)** | external kinematics for `IntegralFamily` (`kinematics(...)` builds it from strings) |

Helpers that only serve other functions (`scalar.V`, `scalar.vlog`, `scalar.vli2`, `ff.Thiele`, `pv.Kin`
and the methods `add_identity`, `clean`, `substitute`) are left out.
