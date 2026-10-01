#import "style.typ": *
#show: setup

#let repo = link("https://github.com/aburousan/feynsage")[github.com/aburousan/feynsage]
#let walkthrough = link("https://github.com/aburousan/feynsage/blob/main/examples/feynsage_walkthrough.ipynb")[`examples/feynsage_walkthrough.ipynb`]
#let reference = link("https://github.com/aburousan/feynsage/blob/main/docs/REFERENCE.md")[`docs/REFERENCE.md`]

// ---------------------------------------------------------------- title page
#page(header: none, numbering: none)[
  #v(1.2cm)
  #align(center)[
    #text(10pt, fill: ink2, tracking: 1.5pt)[A STEP-BY-STEP TUTORIAL WITH CODE AND OUTPUT]
    #v(0.8cm)
    #text(28pt, weight: "semibold", fill: c-blue.darken(15%))[feynsage]
    #v(0.1cm)
    #text(15pt, fill: ink2, style: "italic")[Feynman integrals in SageMath,\ the lecture in code and a real prediction]
    #v(1.0cm)
    #text(12pt)[Kazi Abu Rousan] \
    #text(10pt, fill: ink2)[School of Physical Sciences, NISER, Bhubaneswar] \
    #v(0.35cm)
    #text(11pt)[Source code and install instructions at #repo]
  ]
  #v(1.0cm)
  #dbox("What this tutorial is", c-pink)[
    It teaches the package feynsage by doing with it the calculations of Prof. B. Ananthanarayan's two-part lecture "Feynman integral calculus" (NISER, 2026). These are the tadpole, integration by parts, the bubble and its differential equation, Feynman parameters and the graph polynomials $cal(U)$ and $cal(F)$. Each function is explained in a violet box the first time it is used. The box says what the function takes, what it gives back and what it saves you from doing by hand. The last chapter puts the pieces together for a real prediction (the lifetime of the neutral pion) and compares it with experiment.

    Every grey box marked *In [n]* is real code. It was run, in this order, in one SageMath session (as in a Jupyter notebook) and the box under it marked *Out [n]* is what came out. Nothing in the output boxes was typed by hand. I am also writing up the lecture itself in full (my notes _Feynman integrals at one loop_). They are not finished yet and will be shared once they are complete.
  ]
  #v(0.2cm)
  #text(9pt, fill: ink2)[I wrote feynsage in 2024 for my MSc project. After sir's talk it was refined, checked against the calculations of the lecture and published at #repo under the MIT licence.]
]

#page(header: none)[
  #set text(9.5pt)
  #show outline.entry.where(level: 1): it => { v(6pt, weak: true); strong(it) }
  #outline(depth: 2, indent: 1em)
]

// ======================================================================== 1
= Getting feynsage

== What you need

- *SageMath* 10.7 or newer. It is free. The install script below finds it or installs it.
- *FORM* for the Dirac traces of the last chapter. The script installs it too.
- macOS (Apple silicon or Intel), Linux, or Windows through WSL. SageMath does not run on Windows itself.

== Download and install

The code is at #repo. Either clone it with git or, on the web page, press *Code* and then *Download ZIP* and unpack it. Then, inside the folder:
```bash
git clone https://github.com/aburousan/feynsage.git
cd feynsage
./install.sh            # finds or installs SageMath and FORM, then installs feynsage
./install.sh --test     # the same, then runs the tests
```
If SageMath is already there, `sage -pip install --user --no-deps .` does the last step by hand.

== Check and start

```bash
sage -c "import feynsage; print(feynsage.__file__)"     # where it was installed
sage -n jupyterlab                                       # notebooks, choose the SageMath kernel
```
Every cell of this tutorial can be run again with `sage docs/tutorial/run_cells.sage`.

// ======================================================================== 2
= A short tour of feynsage

feynsage is a Python package that runs inside SageMath, so all of Sage (symbols, exact fractions, calculus, plots) is available next to it. Load it once:
#cell("import")
The first line loads everything that is used most. The second loads a few closed formulas from the module `feynsage.oneloop`. The functions used in this tutorial fall into five groups:

#table(columns: (20%, 50%, 30%),
  table.header([*task*], [*functions*], [*lecture*]),
  [closed formulas], [`tadpole`, `bubble_equal_mass`, `bubble_massless`, `expand_eps`], [tadpole, bubble, the $epsilon$ expansion],
  [families and IBP], [`family`, `.ibp`, `ibp_reduce`, `.draw`, `.UF`], [IBP identities, reduction, completing the square],
  [one-loop functions], [`A0`, `B0`, `C0`, `D0`, `finite_part`, `quick_plot`], [the basis $A, B, C, D$],
  [graphs], [`graph`, `.U`, `.F`, `.spanning_trees`, `.two_forests`, `.plot`, `draw_panels`], [1-trees, 2-trees, $cal(U)$, $cal(F)$],
  [Dirac traces], [`form.dirac_trace`, `form.compute` (traces and contractions in Sage notation, converted to FORM)], [the pion (last chapter)],
)

These are only the functions this tutorial needs. feynsage has many more that the lecture did not reach. Some of them are one-loop tensor integrals in the style of Package-X (`loop`, `PVB`, `PVC`, `PVD`), closed forms of $C_0$ and $D_0$ in dilogarithms (`explicit`, `c0_closed`, `d0_closed`), infrared poles (the module `ir`), the finite-field reducer `reduce_ff`, `symmetries`, `diagram` (ready-made graphs like `"kite"`) and pictures of sectors (`draw_sectors`). All of them are listed with their arguments in #reference in the repository and the notebook #walkthrough uses most of them. `info(f)` explains any one of them.

Two normalisations appear and each is used in its own place.
- The closed formulas and the families use the Euclidean measure of the lecture with the $pi^(D\/2)$ taken out, $integral dif^D k\/pi^(D\/2)$ with propagators $k^2 + m^2$.
- `A0`, `B0`, `C0`, `D0` follow Package-X (Minkowski space, $mu^(2 epsilon) e^(epsilon gamma_E) integral dif^D k\/(i pi^(D\/2))$), so that their numbers can be compared with Package-X and LoopTools.

Every feynsage function explains itself. `info(f)` prints what it computes and what its output means and most functions also take `explain=True`. Sage's own functions are explained with `f?` in a notebook.

// ======================================================================== 3
= Part 1: the tadpole

The lecture started from the simplest loop integral, the tadpole, in Euclidean signature,
$ I_1 = integral (dif^D k)/(k^2 + m^2) $
In $D = 4$ it diverges at large $k$ (ultraviolet). In dimensional regularisation $D = 2 omega$ (Ramond) or $D = 4 - 2 epsilon$ and the angular integral gives the solid angle of $D$ dimensions.

== The solid angle

From the Gaussian integral, $Omega_n = 2 pi^(n\/2) \/ Gamma(n\/2)$. Sage handles such formulas directly:
#func("var('n'),  expr.subs(n=4)", what: "Sage")[
  `var` makes a symbol (here $n$). Formulas built from symbols stay exact. `subs` puts in values. A list in square brackets evaluates the formula several times.
]
For $n = 2, 3, 4$ it must give the circle, the sphere and the 3-sphere:
#cell("solidangle")

== The tadpole with $n$ powers

With $k = m tan t$ the radial integral is a Beta function and
$ I_n = integral (dif^D k)/((k^2 + m^2)^n) = (pi^omega Gamma(n - omega))/(Gamma(n) (m^2)^(n - omega)) $
#func("tadpole(n, m2)")[
  The one-loop tadpole $integral dif^D k\/pi^(D\/2) thin (k^2 + m^2)^(-n)$ in closed form, $Gamma(n - D\/2) (m^2)^(D\/2 - n)\/Gamma(n)$.
  - `n`: the power of the propagator (a number or a symbol)
  - `m2`: the mass squared
  - returns a Sage expression in the symbol `D` (the dimension). The $pi^(D\/2)$ is taken out of the measure, so the lecture's $I_n$ is this times $pi^omega$.
]
#cell("In")
The lecture also showed that $I_n$ follows from $I_1$ by differentiating with respect to $m^2$,
$ I_n (m^2) = (-1)^(n-1)/((n-1)!) (dif/(dif m^2))^(n-1) I_1 (m^2) $
#func("diff(f, x, k),  f.simplify_full()", what: "Sage")[
  `diff` differentiates $k$ times with respect to `x`. `simplify_full` brings an expression to a simple form, so a difference that is zero prints as 0.
]
#cell("derivrel")

== The pole and the $lambda phi^4$ tadpole

Near $D = 4$ the Gamma function has a pole.
#func("expand_eps(expr, order=0, loops=1)")[
  Puts $D = 4 - 2 epsilon$ into an expression with Gamma functions, multiplies by $e^(L epsilon gamma_E)$ ($L$ = `loops`, the usual convention that removes $gamma_E$) and expands in $epsilon$ up to $epsilon^"order"$. Gamma functions at negative integers are rewritten first, so the poles come out correctly.
]
#cell("pole")
This reads $-m^2\/epsilon - m^2 + m^2 log m^2$, a simple pole as the lecture said. In $lambda phi^4$ theory (Ramond, section A.4) the tadpole diagram is $1/2 (-lambda) (mu^2)^(2 - omega) integral dif^(2 omega) ell\/(2 pi)^(2 omega) thin 1\/(ell^2 + m^2)$. Put in the integral, expand around $omega = 2$ with $e = 2 - omega$ and compare with the lecture's result
$ (lambda m^2)/(32 pi^2) {1/(2 - omega) + psi(2) + ln((4 pi mu^2)/(m^2)) + O(2 - omega)} $
#func("expr.series(e, k).truncate()", what: "Sage")[
  The Laurent series of `expr` in `e` up to (not including) $e^k$, as an ordinary expression.
]
#cell("phi4")
The difference is exactly zero. The digamma function at 2 is $psi(2) = 1 - gamma_E$:
#cell("psi2")

// ======================================================================== 4
= Part 1: integration by parts for the tadpole

== Integral families

An *integral family* is a list of propagators $D_1, D_2, dots$ together with all integrals $integral 1\/(D_1^(a_1) D_2^(a_2) dots)$ with integer powers. IBP relates its members. For the tadpole there is one propagator and the members are $T(n) = integral 1\/(k^2 + m^2)^n$.
#func("family(props, kin=None, loops=None, euclidean=False, name=\"F\")")[
  Builds an integral family from strings.
  - `props`: the propagators, each a momentum `"k"` or `("momentum", "mass")` for a massive line
  - `kin`: the scalar products of the external momenta, for example `{"p^2": "s"}`
  - `loops`: the loop momenta (guessed from the names if left out)
  - `name`: the letter used for its members, `"F"` if left out. Here `"T"` for the tadpoles and later `"J"` for the bubbles, as in the lecture. It is only a label.
  - `euclidean=True`: propagators $q^2 + m^2$ as in the lecture (default Minkowski, $q^2 - m^2$)
  - returns a `Family`. A member is written with the family's name and the powers of the propagators. For example `T(2)` is $integral dif^D k\/pi^(D\/2) thin 1\/(k^2 + m^2)^2$ and for a family with two lines `J(2,1)` is $integral 1\/(D_1^2 D_2)$. A power 0 means the line is absent and a negative power is a numerator. (This $F$ or $T$ or $J$ has nothing to do with the hypergeometric function.) `.info()` prints the family, `.ibp(a)` gives the IBP identities, `.UF()` the graph polynomials.
]
#cell("tadfam")

== The IBP identity

This is Gauss's law in $D$ dimensions. The integral of a total derivative vanishes,
$ integral dif^D k thin partial/(partial k^mu) (k^mu)/((k^2 + m^2)^n) = 0 quad => quad 2(omega - n) I_n + 2 n m^2 I_(n+1) = 0 $
#func("fam.ibp(a)")[
  All IBP identities with the *seed* `a` (a tuple of powers, here `(n,)` for $T(n)$). For every loop momentum $ell$ and every vector $v$ (loop or external) it gives the identity from $integral partial\/partial ell^mu thin (v^mu thin "integrand") = 0$, worked out exactly as in the lecture (derivative, then every scalar product written through the propagators). It returns a list of dictionaries `{powers: coefficient}`. For example `{(1,): d - 2, (2,): 2*m^2}` means $(d - 2) T(1) + 2 m^2 T(2) = 0$. `d` is the dimension.
]
#cell("tadibp")
The first line reads $(d - 2) T(1) + 2 m^2 T(2) = 0$, the identity of the lecture with $n = 1$. In general $(d - 2n) T(n) + 2 n m^2 T(n+1) = 0$.

== Solving the identities: reduction

Used again and again, the identities bring every $T(n)$ down to $T(1)$. That is what a reduction program does.
#func("ibp_reduce(fam, targets, method=\"auto\")")[
  Writes every target integral as a combination of *master integrals*, the members that the identities cannot remove.
  - `targets`: the integrals to reduce, by their powers, for example `["T(2)", "T(3)"]` or the same as tuples `[(2,), (3,)]`
  - it chooses which seeds to use, finds the symmetries of the family, writes all IBP identities and solves them (Laporta's method). `method` = `"ff"` (finite fields, fast), `"trimmed"` or `"exact"` (exact rational functions). `"auto"` picks one.
  - returns a `Reduction`. `r[target]` is the dictionary `{master: coefficient}`, `r.masters` the list of masters, `r.info()` a summary. In a notebook it prints as typeset equations.
]
#cell("tadred")
The letters can also be drawn as pictures. `draw` puts the graph of each integral in place of its name. A dot on a line is one extra power of that propagator, so $T(2)$ is the tadpole with one dot.
#func("r.draw(graph, targets=None, rename=None, size=1.5)")[
  Draws the reduction `r` as equations between Feynman diagrams, one row for each target.
  - `graph`: the diagram of the family, made with `graph(...)` (Part 2 explains it). Its lines must be in the same order as the propagators of the family.
  - `targets`: which rows to draw (all of them by default)
  - `rename`: a dictionary to print a variable under another name, for example `{"kk": "k^2"}`
  - a line with power 0 is contracted to a point, so the bubble with one line removed is drawn as a tadpole
]
#cell("taddraw", fig: 58%)
The ratios must be exactly those of the Gamma functions above. The coefficients come back as rational functions of `d`. `SR(str(...))` turns one into an ordinary Sage expression so that it can be compared with `tadpole`:
#cell("tadcheck")

== The difference equation and its solution

Argeri and Mastrolia write the identity for $U(n) = pi^(-D\/2) integral dif^D k\/(k^2 + 1)^n$ as a difference equation, $-(n - D\/2) U(n) + n U(n+1) = 0$. Mathematica's `RSolve` solves it once $U(1)$ is given. The solution is $U(n) = Gamma(n - D\/2) \/ (Gamma(1 - D\/2) Gamma(n)) thin U(1)$. Put it back:
#cell("diffeq")
The lecture also solved it with a *factorial series*, $U(n) = integral_0^1 dif t thin t^(n-1) v(t)$. Putting this into the difference equation and integrating by parts (with $v(1) = 0$) gives a first-order differential equation for $v$,
$ (v'(t))/(v(t)) = (D\/2 - t)/(t(t - 1)) quad => quad v(t) = v_0 thin t^(-D\/2) (1 - t)^(D\/2 - 1) $
and the asymptotic estimate fixes $v_0 = 1\/Gamma(D\/2)$.
#func("desolve(eq, y, ivar=t),  numerical_integral(f, a, b)", what: "Sage")[
  `desolve` solves an ordinary differential equation for the function `y` of `t` (here `_C` is the free constant $v_0$). `numerical_integral` integrates numerically and returns the value and an error estimate. `[0]` takes the value.
]
Sage solves the equation and the integral gives back $U(n)$ (here $D = 2.6$, $n = 3$):
#cell("factorial")

// ======================================================================== 5
= Part 1: the basis of one-loop integrals

Passarino and Veltman showed that every one-loop integral is a combination of four scalar functions, $A$ (tadpole), $B$ (bubble), $C$ (triangle) and $D$ (box). Package-X computes them in Mathematica and FeynCalc uses it through FeynHelpers. FIRE and Kira do the IBP reduction for many loops.
#func("A0(m), B0(s, m1, m2), C0(s1, s12, s2, m0, m1, m2), D0(...)")[
  The four scalar one-loop functions in the Package-X normalisation (Minkowski space). `B0(s, m1, m2)` has propagators $(ell, m_1)$ and $(ell + p, m_2)$ with $p^2 = s$. $A_0$ and $B_0$ come back as formulas with the pole $1\/epsilon$ and the scale $mu$. A finite $C_0$ or $D_0$ stays a symbol until `.n()` asks for its number (more than 30 digits are available).
]
#func("finite_part(expr, mu_value=None),  uv_part(expr)")[
  The $epsilon^0$ part and the coefficient of $1\/epsilon$ of an expression. `mu_value` sets the scale $mu$ to a number.
]
#cell("basis")

// ======================================================================== 6
= Part 1: the bubble

The lecture followed Argeri-Mastrolia and Ramond, with equal masses, external momentum $k$ and loop momentum $p$,
$ J_(n_1 n_2) = integral dif^D p thin 1/(D_1^(n_1) D_2^(n_2)), quad D_1 = p^2 + m^2, quad D_2 = (p - k)^2 + m^2 $
In `family` this is two massive propagators. The name `kk` is just $k^2$. A Sage variable name cannot contain `^`, so $k^2$ is written `kk` (and $k dot p$ would be `kp`):
#cell("bubfam")

== IBP for the bubble

Now there are two vectors to put inside the derivative, the internal momentum $p$ and the external one $k$, so `ibp` gives two identities:
#cell("bubibp")
A dot on a line means one more power of that propagator. `ibp_reduce` solves the identities and gives the two results of the lecture,
$ J_(21) = -(D - 3)/(k^2 + 4 m^2) thin J_(11) + 1/(k^2 + 4m^2) thin T(2), quad quad T(2) = -(D - 2)/(2 m^2) thin T(1) $
where `J(0,1)` is the tadpole $T(1)$ (line 1 absent, power 0):
#cell("bubred")
The coefficients are rational functions in the variables `d`, `kk`, `m`. `.parent().gens()` hands these variables to us, so the comparison with the lecture is exact:
#cell("bubcheck")
The same reduction as pictures. The graph is two lines between the vertices $A$ and $B$, with $k$ coming in at $A$ and going out at $B$. `rename` prints `kk` as $k^2$:
#cell("bubdraw", fig: 100%)

== The differential equation

$J$ depends on $k$ only through $k^2$, so $k_mu partial J\/partial k_mu = 2 k^2 thin partial J\/partial k^2$. Differentiating under the integral,
$ k_mu (partial J)/(partial k_mu) = T(2) - J_(11) - k^2 J_(12) $
The reduction above already contains $J_(12)$ and $T(2)$ in terms of the masters, so the right-hand side becomes $A J + B thin T(1)$. That is the differential equation. Compare it with the lecture,
$ (dif J)/(dif k^2) + 1/2 [1/k^2 - (D - 3)/(k^2 + 4m^2)] J = -(D - 2)/(4 m^2) [1/k^2 - 1/(k^2 + 4m^2)] T(1) $
#cell("bubde")
== Solving the differential equation

Now we solve it, with its boundary condition. Put $k^2 = 4 m^2 x$, so that $dif J\/dif x = 4 m^2 thin dif J\/dif k^2$. `partial_fraction` splits the coefficients into simple poles. The equation has singular points only at $x = 0$ ($k^2 = 0$) and $x = -1$ (the threshold $k^2 = -4m^2$ in this Euclidean metric):
#cell("dex")
#func("expr.subs(kk=...),  expr.subs({dd: D}),  expr.partial_fraction(x)", what: "Sage")[
  `subs` replaces a variable. The keyword form `subs(kk=...)` works when the Python name and the Sage name are the same. `dd` is the Sage variable printed `d`, so it goes in a dictionary. `partial_fraction(x)` writes a rational function of `x` as a sum of simple fractions.
]
*The homogeneous equation.* `desolve` gives $J_0 = C thin x^(-1\/2) (1 + x)^((D - 3)\/2)$:
#cell("dehom")
*The boundary condition.* At $k^2 = 0$ the bubble is $J(D, 0) = integral dif^D p\/(p^2 + m^2)^2 = T(2)$, a finite number. The homogeneous solution behaves as $x^(-1\/2)$ there, so it must not appear. Try a power series $J = sum_n c_n x^n$, which is regular at $x = 0$. Multiply the equation by $x(1 + x)$ to clear the denominators and ask that every power of $x$ cancels. The $x^0$ equation does not contain any free constant. It *fixes* $c_0$ and the value is exactly $T(2)$. So regularity at $x = 0$ and the boundary condition are the same statement. Every higher $c_n$ then follows from the one before:
#cell("deseries")
#func("solve(eqs, vars, solution_dict=True),  expr.coefficient(x, k)", what: "Sage")[
  `coefficient(x, k)` takes the coefficient of $x^k$ in an expanded expression. `solve` solves a list of equations for the listed variables. With `solution_dict=True` it returns a list of dictionaries `{variable: value}`.
]
The ratio of two neighbouring terms of $attach(F, bl: 2, br: 1)(a, b; c; z) = sum_n (a)_n (b)_n\/((c)_n n!) thin z^n$ is $(a + n)(b + n)\/((c + n)(1 + n)) thin z$. The ratios $c_(n+1)\/c_n$ found above are exactly this with $a = 2 - D\/2$, $b = 1$, $c = 3\/2$, $z = -x$. With $c_0 = T(2) = Gamma(2 - D\/2) m^(D - 4)$ the solution is
$ J(D, k^2) = Gamma(2 - D\/2) thin m^(D - 4) thin attach(F, bl: 2, br: 1)(2 - D\/2, 1; 3\/2; -k^2\/4m^2) $
#func("bubble_equal_mass(p2, m2)")[
  The equal-mass bubble $integral dif^D p\/pi^(D\/2) thin 1\/((p^2 + m^2)((p - k)^2 + m^2))$ at $k^2 =$ `p2`, this solution of the differential equation. The kinematic variable sits in the argument and the dimension in the parameters.
]
*A numerical check.* Forget the closed form and integrate the equation with a numerical solver (SciPy's `solve_ivp`), starting next to $x = 0$ from the boundary value. The numbers agree with the hypergeometric function to all ten digits shown ($D = 3.3$, $m = 1$):
#cell("denum")
#func("fast_callable(expr, vars=[x]),  solve_ivp(f, [x0, x1], y0, dense_output=True)", what: "Sage, SciPy")[
  `fast_callable` turns a Sage expression into a fast numerical function. `solve_ivp` integrates $y' = f(x, y)$ from `x0` to `x1` with the starting value `y0`. `num.sol(x)` gives the solution at any `x` in between.
]

== The hypergeometric function in Sage

Sage knows $attach(F, bl: 2, br: 1)$ as `hypergeometric([a, b], [c], z)`. `.n(digits=25)` evaluates it to 25 digits (through mpmath, which Sage contains). Below the same number is found from the series itself and directly from mpmath. The parameters are exact fractions, since a decimal number like `0.35` would limit the precision to 16 digits. For special parameters Maxima can rewrite the function in elementary functions. At $D = 3$ the bubble is an arc tangent:
#cell("hyp")
#func("hypergeometric([a, b], [c], z),  rising_factorial(a, n),  simplify_hypergeometric(algorithm='maxima')", what: "Sage")[
  `hypergeometric` is the generalised hypergeometric function $attach(F, bl: p, br: q)$ with the lists of upper and lower parameters. `rising_factorial(a, n)` is the Pochhammer symbol $(a)_n = a(a + 1) dots (a + n - 1)$. `simplify_hypergeometric` looks for a closed form. Maxima needs to know the sign of $x$, so `assume(x > 0)` comes first.
]

// ======================================================================== 7
= Part 1: the fish diagram and Feynman parameters

The fish diagram of $lambda phi^4$ is the same bubble. The lecture combined the two denominators with the Feynman trick,
$ 1/(A B) = integral_0^1 (dif x)/([x A + (1 - x) B]^2) $
#func("assume(cond),  integrate(f, x, a, b)", what: "Sage")[
  `integrate` does definite integrals exactly. `assume` tells it what it needs to know about the symbols (here that $A, B > 0$).
]
#cell("feyntrick")
After the shift $ell -> ell + p(1 - x)$ the loop integral is a tadpole with mass$""^2$ $m^2 + p^2 x (1 - x)$ and two powers, so `tadpole(2, ...)` gives it. Only the integral over the Feynman parameter $x$ is left. Done numerically it equals the hypergeometric result of the differential equation ($D = 3.3$, $m = 1$, $p^2 = 3$):
#cell("fishfp")

== Ramond's expansion around $D = 4$

Expanding around $omega = 2$ the finite part is $-integral_0^1 dif x ln(1 + x(1-x) p^2\/m^2)$, which Ramond writes as $2 - sqrt(1 + 4m^2\/p^2) thin ln[(sqrt(1 + 4m^2\/p^2) + 1)\/(sqrt(1 + 4m^2\/p^2) - 1)]$:
#cell("ramond")

== Going to Minkowski space: the imaginary part

In Minkowski space $p^2 -> -s$ and above the threshold $s > 4m^2$ the argument of the logarithm turns negative. Then $ln z = ln|z| + i pi$ for $z < 0$ and the bubble gets an imaginary part. By the optical theorem it must, because above threshold the two particles in the loop can be real. `B0` is the same bubble in Minkowski space, so its imaginary part should be $pi beta$ with $beta = sqrt(1 - 4m^2\/s)$ (here $m = 1$):
#cell("continuation")
#func("quick_plot(exprs, (s, a, b), labels=None)")[
  Plots one or more expressions in `s` from `a` to `b`. For an expression with $1\/epsilon$ it plots the finite part. Complex values are drawn as the real part (full line) and the imaginary part (dashed). Returns a matplotlib figure.
]
The imaginary part starts at $s = 4$:
#cell("b0plot", fig: 66%)

// ======================================================================== 8
= Part 2: graph polynomials

The second lecture combined the denominators of any diagram at once and read the result off the graph. Take a diagram with $L$ loops and $N$ propagators $D_j = q_j^2 - m_j^2 + i 0$ (Minkowski space), raised to powers $nu_j$. Its *Feynman integral* is the integral over all loop momenta of the product of propagators. With Feynman parameters it becomes
$ I(nu_1, dots, nu_N) &= integral product_(l=1)^L (dif^D ell_l)/(i pi^(D\/2)) thin 1/(D_1^(nu_1) dots D_N^(nu_N)) \
  &= (-1)^(N_nu) (Gamma(N_nu - L D\/2))/(product_j Gamma(nu_j)) integral_0^infinity product_(j=1)^N dif x_j thin x_j^(nu_j - 1) thin delta(1 - sum x_j) thin (cal(U)^(N_nu - (L+1) D\/2))/(cal(F)^(N_nu - L D\/2)) $
with $N_nu = sum_j nu_j$. This is not yet the scattering amplitude. The amplitude contains such integrals, multiplied by the couplings and by the numerator from the Dirac algebra (the pion in the last chapter shows how they fit together). The sign $(-1)^(N_nu)$ comes from the Minkowski propagators. In Euclidean space, with $dif^D ell\/pi^(D\/2)$ and propagators $ell^2 + m^2$, the same formula holds without it.

The two polynomials come from the graph:
- A 1-tree (spanning tree) is what is left after erasing $L$ lines so that the rest touches every vertex without a loop. Each gives the product of the $x_i$ of the erased lines. $cal(U)$ is their sum, homogeneous of degree $L$.
- A 2-tree is what is left after erasing one more line, so that the graph falls into two pieces. Each gives the product of the erased $x_i$ times the square of the momentum flowing from one piece to the other. Their sum is $V$ and $cal(F) = V + cal(U) sum_i x_i m_i^2$.

#func("graph(edges, legs, kin=None, euclidean=False)")[
  A Feynman diagram as a graph, from strings.
  - `edges`: the lines, `"A-B"` between vertices A and B, or `"A-B:m"` for a line of mass $m$. *The order of the lines is the numbering* $x_1, x_2, dots$ (the lecture's $alpha_1, alpha_2, dots$).
  - `legs`: the momentum *entering* at each vertex, for example `{"A": "p", "B": "-p"}`. They must add up to zero.
  - `kin`: the scalar products of the external momenta. `euclidean=True` for Euclidean signs.
  - returns a `FeynmanGraph` with these methods:
    - `.spanning_trees()`: the 1-trees, as tuples of the *erased* lines (counted from 0)
    - `.two_forests()`: the 2-trees, as pairs (erased lines, vertices of one of the two pieces)
    - `.U()`, `.F0()` ($= V$), `.F()`: the polynomials, read off the trees exactly by the rules above
    - `.plot()`: a picture
    - `.N`, `.V`, `.L`: the numbers of lines, vertices and loops
]
#func("draw_panels(g, erased, titles=None, ncols=5)")[
  Draws the graph `g` several times side by side, each copy with its own erased lines (grey and dashed). Feed it `g.spanning_trees()` to see all the 1-trees.
]

== The bubble with two masses

Two 1-trees (erase the top line or the bottom line) and one 2-tree (erase both lines, then $p$ flows between the two vertices).
#cell("bub2", fig: 74%)
So $cal(U) = alpha_1 + alpha_2$, $V = alpha_1 alpha_2 p^2$ and $cal(F) = V + (m_1^2 alpha_1 + m_2^2 alpha_2) cal(U)$, as in the lecture.

== The box (Smirnov, Fig. 3.6)

The massless box with $p_i^2 = 0$, $s = (p_1 + p_2)^2$ and $t = (p_1 + p_3)^2$, with the lines numbered as on the board. Line 1 is on top (between $p_1$ and $p_3$), 2 on the left (between $p_1$ and $p_2$), 3 on the right (between $p_3$ and $p_4$) and 4 at the bottom (between $p_2$ and $p_4$). feynsage draws the square turned by 45 degrees.
#cell("box", fig: 34%)
All six 2-trees. Only the two that cut the box into a left and a right half, or a top and a bottom half, carry momentum. The other four cut off a single corner, which receives only one $p_i$ with $p_i^2 = 0$, so they do not contribute (the grey pictures of the lecture). The small helper `momentum_in` adds up the momenta entering one piece:
#cell("box2trees", fig: 72%)
#note(title: "The labels matter")[
  With lines 1 and 4 on top and bottom, erasing them separates $p_1, p_2$ from $p_3, p_4$ and gives $s thin alpha_1 alpha_4$. The board had $V = s thin alpha_1 alpha_3 + t thin alpha_2 alpha_4$, which is the same diagram with the lines numbered in order around the box (1 left, 2 top, 3 right, 4 bottom). The safe rule is to look at the picture. $s$ goes with the pair of lines whose removal separates $p_1, p_2$ from $p_3, p_4$.
]

== Smirnov's two-loop example

This is the two-loop propagator diagram of the lecture, with lines 1 and 2 on the left, 3 and 4 on the right and 5 in the middle.
#cell("kite", fig: 34%)
Eight 1-trees, each erasing two lines, give $cal(U)$ of degree 2:
#cell("kitetrees", fig: 100%)
Eight 2-trees carry the momentum $p$, each erasing three lines:
#cell("kite2trees", fig: 100%)
These are the 1-trees and 2-trees drawn in the lecture.

== From $cal(U)$ and $cal(F)$ back to the integral

Put $cal(U)$ and $cal(F)$ into the general formula. For the massless bubble ($N = 2$, $L = 1$) the delta function leaves one integral over $x$. `SR(str(...))` turns the polynomials into ordinary expressions so that $x_1 = x$, $x_2 = 1 - x$ can be put in. At $D = 3.3$ and $p^2 = 1$ the result agrees with the closed form `bubble_massless(1, 1, 1)` (to the precision of the numerical integration):
#cell("fpformula")

== Sir's notebook: completing the square

The Mathematica notebook shown in the lecture finds $cal(U)$ and $cal(F)$ in a second way, without drawing any trees. The formula it uses was not written on the board, so here it is step by step.

*Step 1: one denominator.* With Feynman parameters $x_i$ all the propagators are combined into one sum, $sum_i x_i D_i$. Each $D_i$ is quadratic in the loop momentum, so the sum is too.

*Step 2: sort it by powers of the loop momentum.* For one loop with loop momentum $ell$ every quadratic expression can be written as
$ sum_i x_i D_i = M thin ell^2 - 2 thin Q dot ell + J $
- $M$ collects everything in front of $ell^2$. Every propagator contains $ell^2$ once, so $M = sum_i x_i$.
- $Q$ is a momentum. It collects the terms linear in $ell$ (the factor $-2$ is a convention that makes the next step neat).
- $J$ is the rest, made of external momenta and masses with no $ell$ at all.

*Step 3: complete the square.* Shift $ell = ell' + Q\/M$:
$ M thin ell^2 - 2 Q dot ell + J = M thin ell'^2 + (J - (Q dot Q)\/M) $
Now the integral over $ell'$ is a tadpole, because the loop momentum appears only as $ell'^2$, with a "mass" $(J - Q^2\/M)\/M$. Doing it gives a power of $M$ and a power of $J - Q^2\/M$ and they combine exactly into the general formula above with
$ cal(U) = M = det M, quad quad cal(F) = det M thin (J - Q M^(-1) Q) $
*At more loops* there are several loop momenta $ell_1, ell_2, dots$ The coefficients of $ell_r dot ell_s$ form a matrix $M$ (2 by 2 at two loops) and $Q$ becomes a list of momenta (one per loop). The same formulas hold with $det M$ and the matrix inverse $M^(-1)$. This is why the notebook writes them with a determinant.

As a worked example take the bubble with two masses (Euclidean, $D_1 = k^2 + m_1^2$, $D_2 = (k + p)^2 + m_2^2$). Write the scalar products as symbols (`kk` $= k^2$, `kp` $= k dot p$, `pp` $= p^2$), expand and read off $M$, $Q$ and $J$ as the coefficients.
#cell("square")
This is the same $cal(F)$ as from the 1-trees and 2-trees above.

#func("fam.UF()")[
  $cal(U)$ and $cal(F)$ of a family by exactly this method. It builds $M$, $Q$ and $J$ from the list of propagators and returns $det M$ and $det M (J - Q M^(-1) Q)$, for any number of loops (no graph needed). For a Minkowski family the signs are those of the notebook.
]
#cell("squarecheck")
The notebook's own examples, the triangle, the bubble with two masses and the box (Minkowski signs, as in the notebook), come out as on the screen in the lecture:
#cell("sirnb")

// ======================================================================== 9
= A real prediction: the lifetime of the neutral pion

Now put the tools together for a number that is measured in a laboratory. The neutral pion lives for about $8.4 times 10^(-17)$ s and in $98.8%$ of the cases decays into two photons. It has no charge, so the photons cannot attach to it directly. The decay goes through a *triangle* of charged quarks. The physics is explained below as we go and every step is done with feynsage, with the graph method of Part 2 at its centre.

#physics(title: "The model")[
  The quarks $u$ and $d$ ($Q_u = 2\/3$, $Q_d = -1\/3$) come in $N_c$ colours. They couple to the photon with charge $Q e$ and to the pion through $g thin overline(q) i gamma_5 tau_3 q thin pi^0$ ($tau_3 = +1$ for $u$, $-1$ for $d$), with $g = m\/f_pi$ (the Goldberger-Treiman relation, $f_pi approx 92$ MeV). For each quark there are two triangle diagrams, one for each order of the photons. The photons have $k_1^2 = k_2^2 = 0$ and the pion $(k_1 + k_2)^2 = m_pi^2$.
]

== Step 1: the graph, its trees and its polynomials

Write the triangle as a graph. The pion vertex is P and the photon vertices are A and B. All three lines are quarks of mass $m$. Line 1 is P-A (parameter $x_1$), line 2 is B-P ($x_2$) and line 3 is A-B ($x_3$). The pion brings in $k_1 + k_2$ at P and the photons take out $k_1$ at A and $k_2$ at B (Minkowski space, $k_1^2 = k_2^2 = 0$, $2 k_1 dot k_2 = m_pi^2$):
#cell("piongraph", fig: 30%)
Now the rules of Part 2, drawn. One loop, so a 1-tree erases one line and a 2-tree erases two. Each 2-tree cuts off one vertex. Only the one that cuts off the pion vertex P carries momentum, $(k_1 + k_2)^2 = m_pi^2$. The two that cut off a photon vertex carry $k_1^2 = 0$ or $k_2^2 = 0$ and drop out:
#cell("piontrees", fig: 72%)
So $cal(U) = x_1 + x_2 + x_3$ and $V = x_1 x_2 thin m_pi^2$. With the mass term $cal(F) = -V + cal(U) (x_1 + x_2 + x_3) m^2$ (Minkowski sign), which is the $cal(F)$ printed above. On the simplex $x_1 + x_2 + x_3 = 1$, with $x_1 = x$ and $x_2 = y$, the polynomials become $cal(U) = 1$ and
#cell("piondelta")
This is the $Delta$ of the textbook calculation, read off the trees without completing any square.

== Step 2: the integral from $cal(U)$ and $cal(F)$

Now put $cal(U)$ and $cal(F)$ into the general formula of Part 2. Here $N = 3$ lines, all with power 1, $L = 1$ loop and $D = 4$ (the triangle is finite):
$ integral (dif^4 ell)/(i pi^2) thin 1/(D_1 D_2 D_3) = (-1)^3 thin Gamma(3 - 2) integral_0^infinity dif x_1 dif x_2 dif x_3 thin delta(1 - x_1 - x_2 - x_3) thin (cal(U)^(3 - 4))/(cal(F)^(3 - 2)) = -integral dif^3 x thin delta(dots) thin 1/(cal(U) cal(F)) $
feynsage writes the integrand straight from the polynomials of Step 1:
#cell("pionUF")
That is $-1\/Delta$ on the simplex, so the triangle is $-integral_0^1 dif x integral_0^(1-x) dif y thin 1\/(m^2 - x y m_pi^2) = -I(r)\/m^2$ with
$ I(r) = integral_0^1 dif x integral_0^(1-x) dif y thin 1/(1 - x y r), quad r = (m_pi^2)/(m^2) $
Do the integral numerically from $cal(U)$ and $cal(F)$ (here $m = 1$, $r = 0.5$) and compare with the closed form $I(r) = (2\/r) arcsin^2(sqrt(r)\/2)$:
#cell("pionnum")
The amplitude uses the measure $dif^4 ell\/(2 pi)^4$, which is $i pi^2\/(2 pi)^4$ times the measure of the formula:
#cell("pionnorm")
For $I(r)$ itself, expand $1\/(1 - x y r)$ in powers of $r$, integrate term by term and compare with the closed form.
#cell("pionIr")
For a heavy quark $I -> 1\/2$, the area of the triangle. Even for a light constituent quark ($m = 330$ MeV) $2 I(r)$ is only $1.4%$ above 1.

== Step 3: the Dirac trace

The numerator of the diagram is a trace of Dirac matrices. The pion vertex gives $gamma_5$, each photon vertex a $gamma^mu$ or $gamma^nu$ and each quark propagator its numerator $gamma dot p + m$, where $gamma dot p = gamma^mu p_mu$ (often written with a slash through $p$). With the momenta of the triangle and the loop momentum shifted to $q$ ($ell = q - x k_1 + y k_2$, as in Step 2) the trace is
$ "Tr"[ gamma_5 thin (gamma dot (q - x k_1 - (1 - y) k_2) + m) thin gamma^nu thin (gamma dot (q - x k_1 + y k_2) + m) thin gamma^mu thin (gamma dot (q + (1 - x) k_1 + y k_2) + m) ] $
#func("form.dirac_trace(factors, vectors, dim=4, levi_civita=\"form\")")[
  The trace of a product of Dirac matrices, done by the program FORM. The matrices are given as a list, in the order of the product. Each entry is one matrix:
  #table(columns: (auto, auto), inset: 4pt, stroke: none,
    [`'g5'`], [$gamma_5$],
    [`'mu'` (any name that is not a momentum)], [$gamma^mu$],
    [`'p'` (a momentum)], [$gamma dot p = gamma^mu p_mu$],
    [`'p + m'`], [$gamma dot p + m$ ($m$ times the unit matrix)],
    [`'q - x*k1 + y*k2 + m'`], [$gamma dot (q - x k_1 + y k_2) + m$],
  )
  - `vectors`: the names that are momenta (so that `'q'` means $gamma dot q$ and not an index)
  - `dim=4`: four dimensions, needed with $gamma_5$. A symbol such as `'D'` gives $D$ dimensions
  - `levi_civita="usual"`: write the answer with the usual Levi-Civita symbol $epsilon^(mu nu rho sigma)$ (Peskin and Schroeder, $epsilon^(0123) = -1$), printed $epsilon$. The default `"form"` keeps FORM's own symbol, which is $-i$ times it.
  - returns a Sage expression. `g(mu,nu)` is $g^(mu nu)$ and `dot(p,q)` is $p dot q$. In $epsilon(dots)$ a momentum in a slot means it is contracted there, so $epsilon(k_1, k_2, mu, nu) = epsilon^(rho sigma mu nu) k_(1 rho) k_(2 sigma)$.
]
#cell("piontrace")
So the whole numerator is one term,
$ 4 i m thin epsilon(k_1, k_2, mu, nu) = 4 i m thin epsilon^(mu nu rho sigma) k_(1 rho) k_(2 sigma) $
(moving the two momentum slots past the two index slots is an even permutation). Nothing depends on $q$, $x$ or $y$, so the loop integral of Step 2 is all that is needed. It is a pure scalar triangle. Only the terms with exactly one $m$ survive, because a trace with $gamma_5$ needs an even number, at least four, of other $gamma$ matrices.

== Step 4: the amplitude and why the quark mass drops out

Collect the factors of one quark. They are $N_c$ colours, $(e Q)^2$ from the photon vertices, $g$ from the pion vertex, $4 m$ from the trace (Step 3), $I\/(16 pi^2 m^2)$ from the integral (Step 2) and a factor 2 for the two orders of the photons. Then put in $g = m\/f_pi$ and add $u$ and $d$ (with $tau_3 = plus.minus 1$):
#cell("pionamp")
The $m$ of the coupling cancels the $1\/m$ of the loop. A heavier quark couples more strongly and propagates less. The answer knows only the charges and the number of colours, $A = N_c (Q_u^2 - Q_d^2) e^2\/(4 pi^2 f_pi) = alpha\/(pi f_pi)$ for $N_c = 3$. This is the axial anomaly of Adler, Bell and Jackiw.

== Step 5: the decay rate

Square the amplitude and sum over the photon polarisations. That needs the contraction $epsilon^(mu nu rho sigma) epsilon_(mu nu alpha beta) k_(1 rho) k_(2 sigma) k_1^alpha k_2^beta$. It is written exactly as on paper in Sage notation. `form.compute` turns it into a FORM program, runs it and reads the answer back. No FORM code has to be written (`show_code=True` prints the program it made):
#func("form.compute(expr, vectors, lines=(), rules=None, dim=4, show_code=False)")[
  Traces and index contractions, written in Sage notation and done by FORM.
  - `expr`: a product such as `eps(mu,nu,rho,sigma)*k1(rho)`, with `eps(...)` the Levi-Civita symbol, `g(mu,nu)` the metric, `dot(p,q)` or `p.q` a scalar product, `p(mu)` the component $p^mu$ and any other symbols
  - `lines`: fermion lines to trace, each a list of factors as in `dirac_trace`. For example `[["pp","mu","p","nu"], ["k","mu","kp","nu"]]` is the product of two traces of $e^+ e^- -> mu^+ mu^-$
  - `rules`: scalar products to put in, for example `{"k1.k1": 0, "k1.k2": "s/2"}`
  - `dim`: 4, or a symbol like `'D'` for traces in $D$ dimensions
  - `levi_civita="usual"`: `eps` is the usual $epsilon^(mu nu rho sigma)$, in the input and in the answer (otherwise FORM's own symbol, $-i epsilon$)
  - every repeated index is summed. Returns a Sage expression, with `dot(p,q)` for $p dot q$.
]
#cell("pioneps")
The printed program shows what was done for us. It has FORM's `e_` for $epsilon$ (with the factor $i$ taken care of), the indices declared, `contract` and the two rules $k_1^2 = k_2^2 = 0$. The answer is $epsilon^(mu nu rho sigma) epsilon_(mu nu alpha beta) k_(1 rho) k_(2 sigma) k_1^alpha k_2^beta = 2 (k_1 dot k_2)^2 = m_pi^4\/2$.

#note(title: "Writing FORM yourself")[
  For anything `compute` does not cover, `form.run_form(code)` runs a FORM program written by hand and returns FORM's output. A pedestrian introduction to FORM for particle physics (declarations, `g_` and `trace4`, `e_` and `contract`, `id` and a full Compton-scattering calculation) is on my blog, #link("https://rousan.netlify.app/pages/physics/blogs/form_pedestrian_qft_blog/")[rousan.netlify.app/pages/physics/blogs/form_pedestrian_qft_blog].
]

With the two-body phase space $1\/(8 pi)$, a factor $1\/2$ for the two identical photons and $1\/(2 m_pi)$ for the decaying pion,
$ Gamma = 1/(2 m_pi) dot 1/2 dot 1/(8 pi) dot A^2 (m_pi^4)/2 = (A^2 m_pi^3)/(64 pi) $
Now the numbers, with $alpha = 1\/137.036$, $m_(pi^0) = 134.977$ MeV and $f_pi = 130.2\/sqrt(2)$ MeV (Particle Data Group), against the measurements. The lifetime follows from the width as $tau = ℏ thin "BR"(gamma gamma)\/Gamma$:
#cell("pionrate")

#dbox("Result: the pion lifetime from one triangle", c-blue, breakable: false)[
  #table(columns: 3, inset: 6pt, align: (left, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt + ink2) } else { none },
    [], [*feynsage*], [*experiment*],
    [width $Gamma(pi^0 -> gamma gamma)$], [$7.79$ eV], [$7.80 plus.minus 0.12$ eV (PrimEx-II, 2020): $0.1 sigma$],
    [lifetime from that width], [$8.35 times 10^(-17)$ s], [$8.34 times 10^(-17)$ s],
    [lifetime, world average], [$8.35 times 10^(-17)$ s], [$(8.43 plus.minus 0.13) times 10^(-17)$ s (PDG): $0.6 sigma$],
    [the same with $N_c = 1$], [$0.87$ eV], [ruled out],
  )
]
#physics(title: "Reading the comparison")[
  - The width agrees with the most precise measurement, PrimEx-II, to $0.1 sigma$. The lifetime it implies is the same as ours to the last digit shown.
  - Why 8.35 and not 8.43? The two experimental numbers come from different measurements. The PrimEx-II width alone corresponds to a lifetime of $8.34 times 10^(-17)$ s, which is ours. The PDG lifetime is a world average that also contains older experiments, earlier Primakoff measurements and a direct measurement of the decay length at CERN (1985) that found a longer lifetime. They pull the average to a width of $7.72$ eV, so the average lifetime is a little longer. Our lifetime is $1%$ shorter than this average, which is $0.6$ of its standard deviation. There is no tension.
  - The prediction has its own small uncertainty. Since $Gamma prop 1\/f_pi^2$, so $1%$ in $f_pi$ moves it by $2%$. The corrections left out (the masses of the $u$ and $d$ quarks, mixing with $eta$ and $eta'$) are of a few per cent as well.
  - With one colour the rate would be nine times smaller and the pion would live nine times longer, against every measurement. This decay was one of the early pieces of evidence that quarks come in three colours. The quark mass does not matter (Step 4), which is the content of the anomaly.
]

// ======================================================================== end
= Where to go next

#note(title: "What feynsage is for")[
  feynsage is meant for small and medium problems and for studying. Every step is exact, can be checked and explains itself. For very large reductions Kira and FIRE remain the better choice. I may continue improving feynsage in the future.
]
- #repo: the code and the latest version.
- #walkthrough: a notebook that goes through the whole package.
- My notes of the lecture, _Feynman integrals at one loop_, with every derivation in full. They are not finished yet and will be shared once complete.
- References used in the lecture
  - P. Ramond, _Field Theory: A Modern Primer_ (appendix A)
  - M. Argeri and P. Mastrolia, _Feynman diagrams and differential equations_
  - V. A. Smirnov, _Analytic Tools for Feynman Integrals_ (chapter 3)
  - S. Weinzierl, _Feynman Integrals_
  - T. Rauh's lecture notes
- For the pion
  - S. L. Adler (1969) and J. S. Bell and R. Jackiw (1969)
  - PrimEx-II Collaboration, Science 368 (2020) 506
  - Particle Data Group, Review of Particle Physics
