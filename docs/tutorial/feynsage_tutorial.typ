#import "style.typ": *
#show: setup

#let repo = link("https://github.com/aburousan/feynsage")[github.com/aburousan/feynsage]

// ---------------------------------------------------------------- title page
#page(header: none, numbering: none)[
  #v(1.2cm)
  #align(center)[
    #text(10pt, fill: ink2, tracking: 1.5pt)[A STEP-BY-STEP TUTORIAL WITH CODE AND OUTPUT]
    #v(0.8cm)
    #text(28pt, weight: "semibold", fill: c-blue.darken(15%))[feynsage]
    #v(0.1cm)
    #text(15pt, fill: ink2, style: "italic")[Feynman integrals in SageMath,\ from graphs to numbers]
    #v(1.0cm)
    #text(12pt)[Kazi Abu Rousan] \
    #text(10pt, fill: ink2)[School of Physical Sciences, NISER, Bhubaneswar] \
    #v(0.35cm)
    #text(11pt)[Source code and install: #repo]
  ]
  #v(1.0cm)
  #dbox("How to read this tutorial", c-pink)[
    Every grey box marked *In [n]* is real code. It was run, in this order, in one SageMath session (as in a Jupyter notebook) and the box under it marked *Out [n]* is what came out: printed text, the value of the last line (typeset from Sage's own LaTeX) or a figure. Nothing in the output boxes was typed by hand. To run the same code yourself, see Chapter 1.

    The notation is the one of the lecture note _Feynman integrals at one loop_ (sir's lecture at NISER, 2026): $D = 4 - 2 epsilon$, Feynman parameters $x_i$, the graph polynomials $cal(U)$ and $cal(F)$.
  ]
  #dbox("Beyond the lecture", c-amber)[
    The lecture covered the tadpole and the bubble, IBP, differential equations, Feynman parameters, the graph polynomials $cal(U)$ and $cal(F)$, Mellin-Barnes, sector decomposition, expansion by regions and Kira. feynsage also does things the lecture did not cover. Each of them is explained here from the start, in an amber box like this one, before it is used.
  ]
  #v(0.2cm)
  #text(9pt, fill: ink2)[I wrote feynsage in 2024 for my MSc project. After sir's talk it was refined, checked against the lecture note and against Package-X, LoopTools, FeynCalc, LiteRed and Kira, and published at #repo under the MIT licence.]
]

#page(header: none)[
  #set text(9.5pt)
  #show outline.entry.where(level: 1): it => { v(6pt, weak: true); strong(it) }
  #outline(depth: 2, indent: 1em)
]

// ======================================================================== 1
= Getting feynsage

== What you need

- *SageMath* 10.7 or newer (tested with 10.7 and 10.9). It is free. The install script below finds it or installs it.
- *FORM*, only for the Dirac traces of Chapter 7. The script installs it too.
- *A C compiler* (optional). With one, the finite-field reducer compiles a small kernel the first time it runs. Without one it runs in pure Python with the same answers, about ten times slower.
- macOS (Apple silicon or Intel), Linux, or Windows through WSL. SageMath does not run on Windows itself.

== Download

The code is at #repo. Either clone it with git
```bash
git clone https://github.com/aburousan/feynsage.git
cd feynsage
```
or, without git, open the page in a browser, press *Code* and then *Download ZIP*, unpack it and open a terminal in the unpacked folder `feynsage-main`.

== Install

From inside the folder:
```bash
./install.sh            # finds or installs SageMath and FORM, then installs feynsage
./install.sh --test     # the same, then runs every test (a few minutes)
./install.sh --no-deps  # installs feynsage only, does not touch SageMath or FORM
```
What the script does, step by step:
+ It looks for SageMath: the `sage` command, a conda environment called `sage`, or the macOS app. If there is none it installs it (with Homebrew on macOS, or as the conda environment `sage` from conda-forge on Linux and WSL). This is a large download.
+ It looks for FORM and installs it if it is missing (Homebrew on macOS, the official release binary into `~/.local/bin` on Linux).
+ It installs feynsage into Sage's own Python with `sage -pip`. On the macOS app this is a user install, so the signed app is never changed.
+ It imports the package once to check it and to compile the C kernel.

If you prefer to do it by hand, the third step alone is
```bash
sage -pip install --user --no-deps .
```

== Check and start

```bash
sage -c "import feynsage; print(feynsage.__file__)"     # where it was installed
sage -n jupyterlab                                       # notebooks, choose the SageMath kernel
```
Then open `tutorials/01_first_steps.ipynb`. To update later: `git pull` and `./install.sh --no-deps`.

== The other material

- `tutorials/`: five short Jupyter notebooks with exercises and solutions, the same steps as this tutorial.
- `examples/feynsage_walkthrough.ipynb`: a longer tour, with the neutral pion lifetime and the QED of Chluba's thesis.
- `docs/DETAILS.md`: every module, every check against other programs, speed and limits.
- `docs/tutorial/`: this tutorial. `sage docs/tutorial/run_cells.sage` runs all its cells again and `typst compile docs/tutorial/feynsage_tutorial.typ` makes the PDF.

// ======================================================================== 2
= First steps: the tadpole and the bubble

Load the package. Everything below is in one namespace, like the cells of a notebook.
#cell("import")

== Ask a function what it does

Every public function carries its own explanation. `info(f)` prints it and most functions also take `explain=True`.
#cell("info")

== The tadpole $A_0$

The simplest loop integral has one propagator. `var` makes a symbol in Sage.
#cell("a0")
The $1\/epsilon$ is the ultraviolet pole and $mu$ is the renormalisation scale. Three small functions pull the pieces apart: `uv_part` (the coefficient of $1\/epsilon$), `pole_parts` (the coefficients of $1\/epsilon^2$ and $1\/epsilon$) and `finite_part` (the $epsilon^0$ part, with the value of $mu$ as an optional second argument).
#cell("a0parts")

#beyond("the normalisation of Package-X")[
  The note works in Euclidean space with $[dif ell] = dif^D ell\/pi^(D\/2)$ and propagators $1\/(ell^2 + m^2)$. The one-loop functions of feynsage ($A_0$, $B_0$, $C_0$, $D_0$ and `loop`) follow instead Package-X and LoopTools, the two programs most people compare with:
  $ A_0(m) = mu^(2 epsilon) e^(epsilon gamma_E) integral (dif^D ell_M)/(i pi^(D\/2)) thin 1/(ell_M^2 - m^2 + i 0) $
  in Minkowski space. The two are the same integral seen two ways. The Wick rotation gives $dif^D ell_M = i thin dif^D ell_E$, which cancels the $1\/i$, and $ell_M^2 - m^2 = -(ell_E^2 + m^2)$, which gives one minus sign for each propagator. The factor $e^(epsilon gamma_E)$ removes the $gamma_E$ that every loop produces. So $A_0$ is minus the note's tadpole times $e^(epsilon gamma_E)$:
]
#cell("norm")

== The bubble $B_0$ and its threshold

Two propagators, $(ell, m_1)$ and $(ell + p, m_2)$ with $p^2 = s$. $Lambda$ is the function that Package-X calls `DiscB`, the logarithm that carries the threshold.
#cell("b0")
Below the threshold $s = 4 m^2$ the bubble is real. Above it the two particles in the loop can be real, so it gets an imaginary part. The note (Minkowski chapter) finds $"Im" B_0 = pi beta$ with $beta = sqrt(1 - 4 m^2\/s)$. Check it with $m = mu = 1$:
#cell("b0numbers")
Below threshold the imaginary part is zero (the $10^(-26)$ is rounding), above it is $pi beta$ to every digit. A plot shows the kink at each threshold, $s = (m_1 + m_2)^2 = 4$ and $9$. `quick_plot` draws the real part as a full line and the imaginary part dashed.
#cell("b0plot", fig: 78%)

== An exact identity

At $s = 0$ with equal masses the bubble is the derivative of the tadpole with respect to $m^2$, which gives $A_0(m) = m^2 (B_0(0; m, m) + 1)$. Sage checks it exactly:
#cell("a0b0")

// ======================================================================== 3
= One loop with numerators, the Package-X way

#beyond("tensor integrals and their reduction")[
  In a real amplitude the loop momentum also appears upstairs, from the Dirac algebra, as $ell^mu$, $ell^mu ell^nu$, $ell dot p$ and so on. Passarino and Veltman (1979) showed that every such *tensor integral* is a combination of the scalar integrals $A_0, B_0, C_0, D_0$. The recipe:
  + Write the most general answer allowed by Lorentz symmetry. For the bubble only $p^mu$ and $g^(mu nu)$ exist, so $integral ell^mu\/(D_1 D_2) = p^mu B_1$ and $integral ell^mu ell^nu\/(D_1 D_2) = g^(mu nu) B_(00) + p^mu p^nu B_(11)$.
  + Contract both sides with every vector and with $g_(mu nu)$. On the left, write each $ell dot p$ and $ell^2$ through the propagators: then a propagator cancels and a simpler integral is left.
  + Solve the linear system for the coefficients $B_1$, $B_(00)$, $B_(11)$, ...

  The matrix of this system is built from the scalar products of the external momenta, the *Gram matrix* $G_(i j) = p_i dot p_j$. Solving it divides by $det G$. `loop(numerator, propagators, kin)` does all of this exactly, with $D = 4 - 2 epsilon$ kept until the end (a $D$ in a coefficient times a $1\/epsilon$ pole leaves a finite *rational term*).
]

`loop(numerator, propagator, ...)` takes strings. Each propagator is `[momentum, mass]` and `kin` gives the scalar products. The answer is one line per tensor structure:
#cell("loop1", math-size: 9pt)

== A vector numerator, checked by hand

With $ell^mu$ upstairs the only vector the answer can point along is $p^mu$:
#cell("loopmu", math-size: 8.5pt)
#derivation(title: "The same coefficient by hand")[
  Contract with $p_mu$ and use $2 ell dot p = D_2 - D_1 - (s + m_1^2 - m_2^2)$, where $D_1 = ell^2 - m_1^2$ and $D_2 = (ell + p)^2 - m_2^2$. Each $D_i$ upstairs cancels one propagator and leaves a tadpole:
  $ B_1 = 1/(2 s) [A_0(m_1) - A_0(m_2) - (s + m_1^2 - m_2^2) B_0(s; m_1, m_2)] $
]
`PVB(0, 1, ...)` is Package-X's `PVB[0, 1, ...]`, that is $B_1$. Compare the two at a point above threshold:
#cell("b1check")
Rank two works the same way. The structures are $g^(mu nu)$ and $p^mu p^nu$:
#cell("loopmunu", math-size: 8.5pt)

== Triangles and boxes

#beyond("$C_0$ and $D_0$ as numbers and as formulas")[
  The scalar triangle and box are the next functions after $B_0$. With Feynman parameters (as in the note) the triangle is a two-dimensional integral over the simplex $x_0 + x_1 + x_2 = 1$,
  $ C_0 = -integral dif^3 x thin delta(1 - x_0 - x_1 - x_2) thin 1/(cal(F) - i 0), quad cal(F) = sum_i x_i m_i^2 - sum_(i < j) x_i x_j s_(i j) $
  where $s_(i j)$ is the square of the momentum flowing between lines $i$ and $j$ (in Minkowski space, so $cal(F)$ carries the minus sign). The minus sign in front is the Wick rotation again. 't Hooft, Veltman and later Denner integrated it in closed form: a sum of logarithms and dilogarithms $"Li"_2$ whose arguments are roots of quadratic equations. The hard part is the $+i 0$: it decides on which side of each branch cut every logarithm sits. feynsage keeps every quantity both as an exact expression and as an 80-digit number with a tiny $+i 0$, lets the number pick the branch and compares the two at the end.
]

A finite triangle stays a symbol until a number is asked for. `c0_numeric` integrates the Feynman-parameter form numerically, a second and independent way. `explicit` writes it in logarithms and dilogarithms:
#cell("c0")
To see that nothing is hidden, do the integral of the box above yourself with SciPy at a point where every invariant is spacelike (there $cal(F) > 0$ and no $i 0$ is needed):
#cell("c0hand")
A box, in the LoopTools order of arguments $D_0(s_1, s_2, s_3, s_4; s_(12), s_(23); m_0, dots, m_3)$:
#cell("d0")

== Soft and collinear poles

#beyond("infrared divergences")[
  The ultraviolet pole comes from large loop momenta. With massless particles there is a second kind, the *infrared* divergences: a massless line can carry a very small momentum (*soft*) or move exactly along a light-like external leg (*collinear*). Both make propagators vanish inside the integration region. In dimensional regularisation they also show up as poles in $epsilon$, and when soft and collinear meet they give $1\/epsilon^2$. In a physical cross section they cancel against the emission of real soft photons or gluons (the Bloch-Nordsieck and KLN theorems), which is why they must be kept exactly. feynsage has the six divergent triangles and sixteen divergent boxes of Ellis and Zanderighi (2008) built in.
]

A massless triangle with two light-like legs:
#cell("cir")
Where does the $1\/epsilon^2$ come from? Build the same triangle as a graph. Of its three 2-forests only one carries momentum, because the other two cut off a single light-like leg with $p_i^2 = 0$:
#cell("irgraph")
So $cal(F) = Q^2 x_1 x_3$ and the parametric integrand contains $(x_1 x_3)^(-1-epsilon)$. Each of $integral_0 dif x_1 thin x_1^(-1-epsilon)$ and $integral_0 dif x_3 thin x_3^(-1-epsilon)$ is $-1\/epsilon$ near its end point: one is the collinear region of one leg, the other of the other leg, and the corner where both vanish is the soft region. Together they give the Dirichlet integral $Gamma(1 + epsilon) Gamma(-epsilon)^2\/Gamma(1 - 2 epsilon)$, whose expansion (with $e^(epsilon gamma_E)$) is the last output: $1\/epsilon^2 - pi^2\/12$, as in the note.

== When the Gram determinant vanishes

#beyond("a zero Gram determinant")[
  The tensor reduction divides by $det G$. For a triangle with $p_1^2 = s_1$, $p_2^2 = s_2$ and $(p_1 - p_2)^2 = s_(12)$:
]
#cell("gram")
At the $g - 2$ point of QED ($p_1 = p_2$, photon momentum $q^2 = 0$) the determinant is exactly zero, so the textbook reduction fails there. feynsage then does not divide at all: the lines with the same momentum are grouped and the coefficients come straight from the Feynman-parameter integral.
#cell("gramc")
Package-X gives $C_1 = 1\/(2 m^2)$ and $C_(00) = 1\/4 + (1\/epsilon + log mu^2\/m^2)\/4$, the same.

// ======================================================================== 4
= Diagrams as graphs: $cal(U)$ and $cal(F)$

This chapter follows the second part of the lecture (graph polynomials).

== A diagram in one line

`graph(edges, legs, kin)`:
- `edges`: the lines, `"A-B"` between two vertices or `"A-B:m"` for a line with mass $m$. *The order of the lines is the numbering* $x_1, x_2, dots$ of $cal(U)$ and $cal(F)$.
- `legs`: the momentum *entering* at each vertex. They must add up to zero.
- `kin`: the scalar products of the external momenta. `euclidean=True` for Euclidean signs.

The number of loops is $L = N - V + 1$ (lines minus vertices plus one).
#cell("graphbub", fig: 34%)
In the pictures a massive line is thick and blue and a massless line thin and black. Each line carries its $x_i$.

== The kite: spanning trees

The two-loop kite of the note is ready-made. It has four vertices: `L` and `R` where $p$ enters and leaves, `T` and `B` at the top and the bottom.
#cell("kiteplot", fig: 36%)
A *spanning tree* touches every vertex and has no closed loop. For $L$ loops we must remove $L$ lines. `spanning_trees()` lists the removed lines (counted from 0). In the pictures the removed lines are grey and dashed:
#cell("kitetrees", fig: 100%)
$cal(U)$ is the sum over spanning trees of the product of the $x_i$ of the *removed* lines:
#cell("kiteU")

== The kite: 2-forests

Remove one more line, $L + 1$ in all, and the graph falls into two trees. If the momentum $p$ has to flow from one tree to the other, the 2-forest gives (product of the removed $x_i$) $times p^2$. If $L$ and $R$ sit in the same tree, no momentum flows between the trees and the term is zero.
#cell("kiteforests", fig: 100%)
In Minkowski space $cal(F) = -cal(F)_0 + cal(U) sum x_i m_i^2$, hence the $-s$:
#cell("kiteF", math-size: 8.5pt)

== Three independent checks

feynsage finds the same polynomials in two more ways and compares. The matrix method of the lecture ($cal(U) = det M$ after completing the square) needs momenta on the lines. The third way, Kirchhoff's theorem, was not in the lecture.

#beyond("Kirchhoff's matrix-tree theorem")[
  Kirchhoff (1847, for electric circuits) showed that spanning trees can be counted with a determinant. Give line $i$ the weight $1\/x_i$ and build the *Laplacian* of the graph: on the diagonal, the sum of the weights of the lines at that vertex, and off the diagonal, minus the weight of the line between the two vertices. Strike out one row and the same column. The determinant of what is left is the sum over all spanning trees $T$ of $product_(i in T) 1\/x_i$. Multiplying by $product_i x_i$ turns "lines in the tree" into "lines removed", which is exactly $cal(U)$. For the kite (rows and columns in the order of the vertices printed):
]
#cell("kirchhoff", math-size: 8.5pt)
#cell("kirchhoff2")

#beyond("momenta from the graph")[
  To use the matrix method, or to build an integral family, every line needs a momentum. `family()` does it like this: pick one spanning tree; every line *not* in it (a *chord*) gets its own loop momentum, which flows around the one closed loop that the chord makes with the tree; the external momenta flow through the tree by momentum conservation at each vertex. There are exactly $L$ chords, so there are $L$ loop momenta.
]
With momenta on the lines, the three ways can be compared on graphs of one, two and three loops:
#cell("threeways")

== The two-loop vertex of the note

`vertex2` is the planar two-loop vertex of the sector-decomposition chapter of the note, with the lines numbered as in the note's picture:
#cell("vertex2", fig: 38%)
#cell("vertex2UF")

// ======================================================================== 5
= Sectors

#beyond("sectors")[
  The lecture used sectors in the reduction of the kite without dwelling on them. They are the bridge between the graph and integration by parts, so this chapter goes slowly: what an index vector is, what a sector is, why a missing line is drawn shrunk to a point, which sectors are zero, and why all this matters for the master integrals.
]

== Index vectors

Take the kite as a *family* of integrals. `family()` gives every line a momentum (as explained above):
#cell("kitefam")
Every member of the family is fixed by five integers, one per line:
$ F(a_1, dots, a_5) = integral dif^D ell_1 dif^D ell_2 thin 1/(D_1^(a_1) D_2^(a_2) D_3^(a_3) D_4^(a_4) D_5^(a_5)) $
- $a_i = 1$: line $i$ is an ordinary propagator
- $a_i = 2, 3, dots$: the same propagator raised to a higher power (a *dot* on the line)
- $a_i = 0$: the propagator is *absent*
- $a_i < 0$: $D_i$ sits upstairs as a *numerator*

== What a sector is

The *sector* of $F(a_1, dots, a_5)$ is the pattern of lines that are really there: a 1 where $a_i > 0$ and a 0 where $a_i lt.eq 0$. Dots and numerators do not change the sector:
#cell("sectorof")
So the sector tells us *which diagram* we are looking at. The dots and numerators only change what sits on its lines. A family with $t$ lines has $2^t$ sectors, here $2^5 = 32$.

== Why a missing line is drawn shrunk to a point

When $a_5 = 0$ the propagator $1\/D_5$ is gone, but the integral is still over both loop momenta. Look at the momenta above: without $D_5$ the remaining four depend on $ell_1$ (lines 1, 2) or on $ell_2$ (lines 3, 4), never on both. So the integral splits into two one-loop bubbles that touch at one point. In the picture this means the two ends of line 5 (the vertices `T` and `B`) become *one vertex*: line 5 is shrunk to a point (contracted). It is *not* erased, since an erased line would leave a one-loop graph with four lines.
#cell("contract", fig: 62%)
The same holds for the graph polynomials. The parametric formula of a sector is the one of the full family with $x_j = 0$ for every missing line $j$ (an index $a_j = 0$ puts $x_j^(a_j - 1)\/Gamma(a_j)$, which becomes $delta(x_j)$). Sage checks that this is exactly $cal(U)$ and $cal(F)$ of the shrunk graph:
#cell("contractUF")
#physics(title: "In one line")[
  A missing line costs nothing to cross. Its two ends act as one point and the diagram of the sector is the original one with that line shrunk.
]

== All sectors with four and three lines

`draw_sectors` draws each sector this way. The title is the sector and *zero* marks the sectors whose integrals all vanish (explained next).
#cell("allsectors", fig: 100%)
How to read them:
- The five four-line sectors are all non-zero. Shrinking line 5 gives the two bubbles in a row, as above. Shrinking any other line gives a triangle with a bubble on one side.
- Of the ten three-line sectors only two survive. Both are *sunsets*: three lines between the vertex where $p$ enters and the vertex where it leaves.
- In every other three-line sector a line closes on itself at one vertex (a little loop), or both legs end up at the same vertex.

== Zero sectors

Take the sector $(1,1,1,0,0)$ and also $(1,0,1,0,1)$:
#cell("zerosector", fig: 60%)
In $(1,1,1,0,0)$ line 3 has become a loop that starts and ends at one vertex. Its loop momentum $k$ appears in this line only, so the integral contains the factor
$ integral dif^D k thin 1/(k^2) $
This has no scale at all: no mass and no external momentum. Under $k -> lambda k$ it changes by $lambda^(D - 2)$, so it can only be zero in dimensional regularisation (the note shows this for the massless tadpole). One zero factor makes the whole sector zero. In $(1,0,1,0,1)$ both legs end up on the same vertex, where $p$ comes in and goes out again. Their momenta cancel there, so the picture has no legs left: no momentum flows through the sunset, which is again a massless vacuum integral and zero.

#beyond("Lee's criterion for zero sectors")[
  The computer cannot look at pictures, so it uses a test of R. N. Lee (2013) on $G = cal(U) + cal(F)$ of the sector.

  *The test.* The sector is zero if there are numbers $k_i$ such that the rescaling $x_i -> lambda^(k_i) x_i$ multiplies *every* monomial of $G$ by the same factor $lambda$.

  *Why it works.* In the Lee-Pomeransky form of the parametric integral the integrand is $G^(-D\/2)$ times powers $x_i^(a_i - 1)$. Change variables $x_i -> lambda^(k_i) x_i$. The integral cannot change (it is only a change of variables), but the integrand picks up a factor $lambda^(c)$ with $c = sum_i k_i a_i - D\/2$, which is not zero for general $D$. A number equal to $lambda^c$ times itself for every $lambda$ is zero.

  *In practice* each monomial $x^e$ of $G$ gives one linear equation $sum_i k_i e_i = 1$, and the sector is zero exactly when these equations have a solution. That is a few lines of linear algebra.
]
For $(1,1,1,0,0)$, $G = x_1 x_3 + x_2 x_3 - s x_1 x_2 x_3$ and $x_3 -> lambda x_3$ does it ($k = (0, 0, 1)$): that is the lonely line 3. For the sunset no such rescaling exists:
#cell("leecrit", math-size: 8.5pt)
Counted by the number of lines:
#cell("count")
Of the 31 sectors with at least one line only 8 are non-zero.

== Why sectors matter

#result(title: "Three facts that make sectors useful")[
  + *IBP never creates a line.* An IBP identity can add a dot or a numerator, and a numerator can cancel a propagator, but no identity puts back a propagator that was absent. So an integral of one sector is always written through integrals of the *same or smaller* sectors. This is why Laporta's algorithm starts with the simplest sectors and works upwards.
  + *Zero sectors drop out at once.* Here 23 of the 31 sectors never need any work.
  + *Every master integral belongs to one sector.* For the kite the masters are one integral in the sunset sector and one in the bubble-bubble sector:
]
#cell("kitemasters", fig: 52%)
#cell("kitered", math-size: 9pt)
The kite itself, the top sector $(1,1,1,1,1)$, has *no* master of its own: everything in it reduces to the two simpler sectors.

== The same sector under another name

The kite has two non-zero sunset sectors, $(0,1,1,0,1)$ and $(1,0,0,1,1)$. Drawn, they are the same diagram: three massless lines between the in-vertex and the out-vertex. So their integrals are equal after renaming the lines and the reducer keeps only one of them. `canon` shows the representative it picks for a few integrals:
#cell("sectorsym", fig: 56%)

#beyond("sector symmetries")[
  *Pak's criterion.* An integral *without numerators* depends only on its powers $a_i$ and on $G = cal(U) + cal(F)$ of its sector (that is all the parametric formula contains). So if renaming the $x_i$ turns the $G$ of one sector into the $G$ of another, the two sectors hold the same integrals, line for line. Here are the two sunsets:
]
#cell("pak")
#note(title: "How the computer compares two polynomials")[
  Renaming $x_2, x_3, x_5 -> x_1, x_4, x_5$ turns one $G$ into the other. To find such renamings for every pair of sectors at once, feynsage turns each $G$ into a small graph (a vertex for each $x_i$, a vertex for each monomial coloured by its coefficient, an edge whenever $x_i$ appears in a monomial) and asks Sage for its *canonical labelling*: a fixed way of numbering the vertices such that two graphs get identical labels exactly when they are the same up to renaming. Sectors with the same canonical form are the same; one of them is chosen as the representative and every integral is moved there. To keep this fast for large families, a cheap fingerprint (number of monomials, their coefficients and degrees) is compared first and only matching sectors are labelled.
]
For the kite the same renaming also follows from a symmetry of the whole diagram (the left-right reflection). Sometimes two sectors match *only* as sectors, with no symmetry of the full family behind them. The two-loop sunset with three equal masses (with the two numerator lines $k + p$ and $ell$ that IBP needs) shows the difference:
#cell("sunset63")
With sector symmetries the three tadpole-product sectors become one and the sunset sector keeps two masters: 3 in all, the same as Kira 3.1 finds (`tests/test_kira_symmetries.sage` compares every coefficient).

#beyond("integrals with numerators")[
  With a numerator Pak's criterion is not enough: renaming the lines does not tell what happens to $D_j$ upstairs. Then feynsage looks for a *shift of the loop momenta* $ell -> A ell + B p$ (with $det A = plus.minus 1$, so the measure does not change, and the external momenta allowed to swap or change sign when their scalar products stay the same) that carries every propagator of the sector into the corresponding propagator of the other. Applied to the numerator, the shift turns $D_j$ into a combination of propagators, which is written out. The result is one more linear equation between integrals, added to the IBP system. For example, for the equal-mass sunset, a symmetry inside the sunset sector gives
]
#cell("numrel")
These equations are added only when a master with numerators could need them; otherwise the reduction is left as it is.

// ======================================================================== 6
= Integration by parts and master integrals

This chapter follows the IBP part of the lecture and the course notebooks.

== The bubble family

The equal-mass bubble of the note, Euclidean, $m = 1$. A massive line is written `("momentum", "mass")`.
#cell("bubfam")
`ibp(a)` gives the identities for the seed $F(a)$, one for each vector ($ell$ and $p$) inside the derivative. Each line reads $sum c thin F = 0$ with $d$ the dimension:
#cell("bubibp")
With $F(1,2) = F(2,1)$ the first one is the identity of the note, $0 = (D - 3) J(1,1) + (p^2 + 4 m^2) J(2,1) - T(2)$.

== Reduction to master integrals

`ibp_reduce(family, targets)` picks the seeds, finds the symmetries and solves the system:
#cell("bubred", math-size: 8.5pt)
The note has $J(2,1) = -(D-3)\/(p^2 + 4m^2) thin J(1,1) + T(2)\/(p^2 + 4m^2)$ and $T(2) = -(D-2)\/(2m^2) thin T(1)$. With $m = 1$ both coefficients agree exactly:
#cell("bubcheck")

== A two-loop number from two masters

The kite in the labelling of the note, $q^2 = 1$. Both masters are products of Gamma functions: the bubble times bubble is $G(1,1)^2$ and the sunset is $G(1,1) G(1, 2 - D\/2)$, with $G(a,b)$ the massless one-loop bubble. Put them into the reduction and expand:
#cell("kitevalue")
All the poles and all the rational numbers cancel and $6 zeta(3)$ is left, the value of the note found there with Feynman parameters and Cheng-Wu.

== Three reducers, one answer

#beyond("finite fields")[
  An IBP system for a real problem has hundreds of thousands of equations whose coefficients are rational functions of $d$ and the invariants. Gaussian elimination with such coefficients swells: intermediate expressions become huge although the final answer is small. Modern programs (FIRE, Kira, FiniteFlow) avoid this:
  + Put numbers for $d$ and the invariants, and do the arithmetic *modulo a large prime* $p$. Every number then fits in one machine word and nothing swells.
  + The first such run shows which equations the targets really need. All the others are thrown away (*trimming*).
  + Repeat at many values of $d$ and rebuild each coefficient as a rational function of $d$ from these values (Thiele interpolation, a continued fraction). With two primes, rebuild the rational numbers inside it (*rational reconstruction*) and check with a third.

  feynsage does this in `method="ff"`, with a small C kernel for the first run. `method="trimmed"` uses the first run only to find the needed equations and then solves those with exact rational functions. `method="exact"` solves everything exactly, which is simple but slow. All three must agree:
]
#cell("methods")
With `verbose=True` the finite-field reducer says how big the system was:
#cell("ffverbose")

== Seeds

`ibp_reduce` also chose the *seeds*: the integrals $F(a)$ for which the IBP identities were written. It takes every non-zero sector, adds up to `rmax` extra dots and `smax` numerator powers (enough for the targets) and orders the integrals Laporta's way: fewer lines first, then fewer dots and numerators. Each identity is solved for its most complicated integral. Whatever cannot be eliminated is a master integral. The note's Laporta chapter shows this by hand.

== Counting masters without reducing

#beyond("the Lee-Pomeransky count")[
  Lee and Pomeransky (2013) found that the number of master integrals can be read off $G = cal(U) + cal(F)$ without any reduction: in each sector, set the missing $x_j$ to zero and count the points where all derivatives of $G$ vanish with every remaining $x_i != 0$. That number is the number of masters of the sector.

  *By hand for the bubble* ($m = 1$): $G = x_1 + x_2 + x_1 x_2 p^2 + (x_1 + x_2)^2$. The two derivatives differ by $(x_2 - x_1) p^2$, so $x_1 = x_2 = x$ and $1 + x(p^2 + 4) = 0$: one point. In the tadpole sector $G = x_1 + x_1^2$ has one point too.

  *Two cautions* (the master chapter of the note has the details): the count does not know about symmetries, so equal sectors must be counted once, and it needs the points to be isolated, which fails in some massless sectors.
]
A Gröbner basis counts the points. The extra variable $t$ with $t product x_i = 1$ throws away the points where some $x_i = 0$:
#cell("lp")
Top sector 1, each tadpole sector 1. The two tadpole sectors are the same integral, so $1 + 1 = 2$ masters, as IBP found.

// ======================================================================== 7
= Dirac traces with FORM

#beyond("FORM")[
  FORM (J. Vermaseren) is the program in which most large perturbative calculations are done: it handles expressions with millions of terms. feynsage writes a short FORM program for each trace, runs it and reads the result back into Sage. Two conventions to know:
  - `'g5'` is $gamma_5$. With $gamma_5$ the trace is done in 4 dimensions (FORM's `trace4`). In $D$ dimensions $gamma_5$ has no unique definition and one must choose a scheme, so feynsage does not mix $gamma_5$ with $D$.
  - FORM's Levi-Civita symbol `e_` is $-i$ times the usual $epsilon^(mu nu rho sigma)$. feynsage prints it as `eps` and keeps FORM's convention.
]

The metric comes back as `g(mu, nu)`:
#cell("tr1", math-size: 9pt)
A repeated index is summed ($gamma^mu gamma_mu = D$ in $D$ dimensions). A factor like `'p + m'` means $slash(p) + m$:
#cell("tr2")

== A whole cross section: $e^+ e^- -> mu^+ mu^-$

For massless fermions summed over spins $sum |cal(M)|^2 = (e^4\/s^2) "Tr"[slash(p)' gamma^mu slash(p) gamma^nu] thin "Tr"[slash(k) gamma_mu slash(k)' gamma_nu]$. Two traces share $mu, nu$, so we write the FORM program ourselves and run it with `form.run_form`. `g_(1, ...)` and `g_(2, ...)` are the two fermion lines. The `id` lines put in the Mandelstam variables, $p dot k = p' dot k' = -t\/2$ and $p dot k' = p' dot k = -u\/2$:
#cell("ee")
Average over the four spin states, go to the centre-of-mass frame ($t = -s(1 - cos theta)\/2$, $u = -s(1 + cos theta)\/2$) and use $dif sigma\/dif Omega = overline(|cal(M)|^2)\/(64 pi^2 s)$:
#cell("eexs")
Since $2 - sin^2 theta = 1 + cos^2 theta$ this is $dif sigma\/dif Omega = alpha^2 (1 + cos^2 theta)\/(4 s)$ and the total cross section is $4 pi alpha^2\/(3 s)$, as in Peskin and Schroeder (eq. 5.12).

// ======================================================================== end
= Where to go next

#note(title: "What feynsage is for")[
  feynsage is meant for small and medium problems and for studying: every step is exact, can be checked and explains itself. On such problems it runs in the same range as the established tools and sometimes a little faster. For very large reductions (millions of equations, several scales) Kira and FIRE remain the better choice. I may continue improving feynsage in the future.
]

- #repo: the code, the issues page and the latest version.
- `tutorials/`: the same steps as five Jupyter notebooks, with exercises and solutions.
- `examples/feynsage_walkthrough.ipynb`: the longer tour, with the neutral pion lifetime and the QED of Chluba's thesis.
- `docs/DETAILS.md`: every module, every check against other programs and the speed numbers.
- The lecture note _Feynman integrals at one loop_: the physics behind every step.
