#import "theme.typ": *
#import "diagrams.typ": *
#show: setup

// ---------------------------------------------------------------- title
#page(footer: none)[
  #v(1.6cm)
  #text(34pt, weight: "semibold", fill: c-blue.darken(25%))[feynsage]
  #v(-0.4em)
  #text(19pt, fill: ink2)[Feynman Integrals in SageMath, from Graphs to Numbers]
  #v(-0.2em)
  #line(length: 70%, stroke: (paint: c-pink.lighten(20%), thickness: 1pt, dash: "dashed"))
  #v(0.6em)
  #text(14pt)[Kazi Abu Rousan] \
  #text(12pt, fill: ink2)[Built on the notes of Prof. B. Ananthanarayan's lecture "Feynman integrals at one loop", NISER, 29 September 2026]
  #v(1.2em)
  #text(12pt, fill: ink2)[github.com/aburousan/feynsage #h(0.6em) · #h(0.6em) MIT licence #h(0.6em) · #h(0.6em) every feynsage result in this talk was computed by the package]
]

// ---------------------------------------------------------------- why
#slide[Why another package?][
  #grid(columns: (1fr, 1fr), gutter: 0.9cm,
  [
    Today a one-loop or multi-loop calculation usually means several programs:
    - Package-X, FeynCalc and LiteRed need Mathematica; FIRE uses it for set-up around a C++ reducer
    - LoopTools gives numbers, not formulas
    - Kira and FIRE reduce, but do not draw or explain

    We wanted *one* free tool where every step is exact and can be checked.
  ],
  [
    #result(title: "feynsage")[
      - written natively for SageMath (free, Python)
      - exact algebra in $QQ(D, s, t, m, dots)$ and graph theory come built in
      - one loop like Package-X, many loops like LiteRed / Kira
      - FORM for the Dirac algebra
      - every function explains its output: `explain=True`
    ]
  ])
]

// ---------------------------------------------------------------- map
#slide[What is inside][
  #set text(13.5pt)
  #table(columns: (24%, 42%, 34%),
    table.header([*Step*], [*feynsage*], [*Plays the role of*]),
    [diagram], [`graph("A-B, A-B:m", legs, kin)`, `diagram("kite")`, pictures], [pen and paper],
    [#U, #F], [spanning trees, 2-forests, Kirchhoff, $det M$], [the speaker's notebook, FeynCalc],
    [family, sectors], [`family([...])`, zero sectors, symmetries], [LiteRed],
    [IBP], [Laporta; finite fields + rational reconstruction], [FIRE, Kira],
    [one loop], [`loop(...)` → $A_0, B_0, C_0, D_0$], [Package-X, FeynCalc],
    [values], [closed forms of $C_0$, $D_0$; IR poles], [Package-X, LoopTools],
    [Dirac traces], [`form.dirac_trace(...)`, $gamma_5$ included], [FORM, FeynCalc],
  )
]

// ---------------------------------------------------------------- start
#slide[Three lines to start][
  ```python
  from feynsage import *
  loop("l^mu l^nu", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})     # one loop, exact
  ibp_reduce(family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1}),
             ["F(2,2,1,2,2)"])                                          # two loops, IBP
  diagram("kite").plot()                                                # a picture
  ```
  #v(0.3em)
  Install: `git clone ... && ./install.sh` #h(0.4em) (macOS, Linux, WSL; it finds or installs Sage and FORM).

  #note(title: "Notation, as in the lecture note")[
    $D = 4 - 2 epsilon$, loop momentum $ell$, Feynman parameters $x_i$, graph polynomials #U and $#F = #F _0 + #U sum_i x_i m_i^2$.
    One-loop functions in the Package-X normalisation $mu^(2 epsilon) e^(epsilon gamma_E) integral dif^D ell \/ (i pi^(D\/2))$.
  ]
]

// ---------------------------------------------------------------- UF: the lecture idea
#slide[#U and #F: the idea of the lecture][
  #grid(columns: (1.15fr, 1fr), gutter: 0.7cm,
  [
    The speaker's notebook (`uf-new-short.nb`) completes the square, loop by loop:
    $ sum_i x_i D_i = sum_(r,s) ell_r dot ell_s M_(r s) - 2 sum_r ell_r dot Q_r + J $
    $ #U = det M, quad quad #F = det M thin (J - Q^T M^(-1) Q) $
    and gives four examples. feynsage reproduces all four, *two ways*:
    #table(columns: 3,
      table.header([*diagram*], [*from propagators*], [*from the graph*]),
      [massless triangle], ok, ok, [bubble, two masses], ok, ok,
      [massless box], ok, ok, [two-loop sunset], ok, ok)
  ],
  [
    #physics(title: "Read off the graph")[
      - #U: cut $L$ lines so that a spanning tree is left, multiply their $x_i$
      - #F#sub[0]: cut $L + 1$ lines so that two trees are left, times the momentum$""^2$ flowing between them
      - #U again: $det$ of the Kirchhoff matrix (matrix-tree theorem)
    ]
    #align(center)[#image("out/kite_trees.svg", height: 4.1cm)]
    #v(-0.4em)
    #align(center)[#text(10pt, fill: ink2)[the 8 spanning trees of the kite, drawn by feynsage]]
  ])
]

// ---------------------------------------------------------------- UF: extended
#let ufrows = json("out/uf_rows.json")
#let fcrows = json("out/uf_feyncalc.json")
#slide[#U and #F: the same idea, pushed further][
  The notebook stops at two loops. Here are bigger diagrams, each done *three independent ways* inside feynsage, and then by FeynCalc's `FCFeynmanPrepare`:
  #v(0.2em)
  #table(columns: (23%, 6%, 6%, 9%, 9%, 11%, 12%, 12%, 12%),
    table.header([*diagram*], [$L$], [$N$], [*trees*], [*terms of #F*], [*trees = $det M$*], [*Kirchhoff = trees*], [*= FeynCalc*], [*time (trees)*]),
    ..ufrows.map(r => (
      [#r.name], [#r.L], [#r.N], [#r.trees], [#r.Fterms],
      if r.matrix { ok } else { bad }, if r.kirchhoff { ok } else { bad },
      if fcrows.at(r.name, default: false) { ok } else { bad },
      [#if r.t_graph < 0.01 [< 10 ms] else [#calc.round(r.t_graph * 1000) ms]],
    )).flatten()
  )
  #v(0.2em)
  #text(11pt, fill: ink2)[The "tennis court" we first tried had a 2-valent vertex (two lines with the same momentum): FeynCalc merges such lines into one doubled propagator, feynsage keeps two parameters. Both are right; the triple box above is the clean test.]
  #grid(columns: (auto, 1fr), gutter: 0.4em, align: horizon,
    [kite: $#U =$], mt("out/kite_U.tex", size: 11pt),
    [$#F =$], mt("out/kite_F.tex", size: 11pt))
]

// ---------------------------------------------------------------- IBP
#slide[Two loops: IBP reduction of the kite][
  #grid(columns: (1fr, 0.42fr), gutter: 0.6cm,
  [
    ```python
    kite = family(["l1", "l1 + q", "l1 + l2", "l1 + l2 + q", "l2"], kin={"q^2": 1})
    ibp_reduce(kite, ["F(2,2,1,2,2)"])
    ```
    #mt("out/kite_red.tex", size: 11pt)
    Two master integrals: the product of two bubbles $F(1,1,1,1,0)$ and the sunset $F(0,1,1,0,1)$. The coefficients are *identical* to LiteRed and Kira 3.1.
  ],
  [
    #image("out/kite.svg", width: 100%)
  ])
  #table(columns: (34%, 18%, 16%, 16%, 16%),
    table.header([*kite, seeds up to*], [*equations*], [*exact Laporta*], [*finite fields*], [*Kira 3.1*]),
    [4 dots, 2 numerators], [9 368], [7.3 s], [0.2 s], [],
    [9 dots, 4 numerators], [157 843], [too slow], [2.7 s], [6.7 s#super[\*]],
  )
  #text(10pt, fill: ink2)[\*Kira on a server with one thread, feynsage on a laptop: a rough guide only.]
]

// ---------------------------------------------------------------- one loop
#slide[One loop, the Package-X way][
  The rank-2 bubble with all symbols, #raw(read("out/bubble_time.txt").trim()) s:
  ```python
  loop("l^mu l^nu", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})
  ```
  #grid(columns: (auto, 1fr), gutter: 0.4em, align: horizon,
    [$g^(mu nu)$:], mt("out/bubble_g.tex", size: 11.5pt),
    [$p^mu p^nu$:], mt("out/bubble_pp.tex", size: 11.5pt))
  #text(10.5pt, fill: ink2)[$Lambda(s; m, m)$ is Package-X's `DiscB`.]
  #grid(columns: (1fr, 1fr), gutter: 0.6cm,
  [
    - projection on a symmetric basis, Gram matrix inverted *exactly* (adjugate)
    - rational terms from $D = 4 - 2 epsilon$, kept to $cal(O)(epsilon^2)$ when IR poles appear
    - `DiscB`, `LogM` carry the $+i epsilon.alt$ branch, as in Package-X
  ],
  [
    - `PVB`, `PVC`, `PVD`: Package-X's coefficient functions
    - Euclidean or Minkowski: `euclidean=True`
    - numerators from FORM traces go straight in
  ])
]

// ---------------------------------------------------------------- closed forms
#slide[$C_0$ and $D_0$ in closed form][
  `explicit(C0(0, 0, 5, 1, 1, 1))`, exact, with the branches fixed by the $+i epsilon.alt$:
  #mt("out/c0_closed.tex", size: 12pt)
  #grid(columns: (1.1fr, 1fr), gutter: 0.6cm,
  [
    `.n(digits=30)` gives #raw(read("out/c0_num.txt").trim().replace("*I", " i"))

    $D_0$: Denner's 16 dilogarithms (no massless line), Denner–Dittmaier for one to four massless lines. `explicit(D0(1,2,3,4,-5,-6,1,1,1,1))` is exact, #read("out/d0_len.txt").trim() characters long, value #raw(read("out/d0_num.txt").trim()).
  ],
  [
    #result(title: "How every branch is fixed")[
      Each quantity is carried twice: as an 80-digit number with a tiny $+i epsilon.alt$ and as an exact expression. Where an argument lies on a cut, the number decides the side. At the end the two are compared at two precisions.
    ]
  ])
]

#slide[Through the thresholds][
  ```python
  quick_plot([B0(s, 1, 1), C0(0, 0, s, 1, 1, 1)], (s, -4, 12), parts="both")
  ```
  #align(center)[#image("out/threshold.svg", height: 72%)]
]

// ---------------------------------------------------------------- IR
#slide[Soft and collinear poles][
  All 6 triangles and 16 boxes of Ellis and Zanderighi, any order of the lines, in the Package-X normalisation:
  #grid(columns: (1fr, 1fr), gutter: 0.5cm,
  [
    `C0(0, 5, 0, 0, 0, 0)` #h(0.3em) (massless, $s = 5$, $mu = 1$):
    #mt("out/ir_tri.tex", size: 11pt)
  ],
  [
    `D0(0,0,0,0, -2,-3, 0,0,0,0)` #h(0.3em) (massless box, $mu = 1$):
    #mt("out/ir_box.tex", size: 11pt)
  ])
  QED vertex $C_0(m^2, q^2, m^2; 0, m, m)$, photon massless, electrons on shell ($m = 1$, $q^2 = -3$): #h(0.3em) $1\/epsilon$ coefficient $0.3419036239$, finite part $-0.1252442693$.
  #note(title: "Checked")[
    118 points against Package-X (all three orders in $epsilon$, double precision). Boxes 15 and 16, whose finite parts Package-X cannot evaluate, against a photon mass $lambda$: $log lambda^2 -> 1\/epsilon + log mu^2$.
  ]
]

// ---------------------------------------------------------------- Gram
#slide[When the Gram determinant vanishes][
  The $g - 2$ vertex at $q^2 = 0$: $p_1 = p_2$, so the usual reduction divides by zero. feynsage takes the coefficients straight from Feynman parameters (lines with identical kinematics are grouped and integrated out exactly):
  #grid(columns: (1fr, 1fr), gutter: 0.6cm,
  [
    ```python
    PVC(0, 1, 0, m^2, 0, m^2, 0, m, m)
    PVC(1, 0, 0, m^2, 0, m^2, 0, m, m)
    ```
    #grid(columns: (auto, 1fr), gutter: 0.4em, align: horizon,
      [$C_1 =$], mt("out/gram_c1.tex"), [$C_(00) =$], mt("out/gram_c00.tex"))
  ],
  [
    #result(title: "Package-X")[
      $C_1 = 1\/(2 m^2)$, #h(0.5em) $C_(00) = 1/4 + 1/4 (1/epsilon + log mu^2/m^2)$ #h(0.3em) #ok
    ]
    Also all invariants zero (divided differences, any masses) and a massive photon: 11 more points agree with Package-X.
  ])
]

// ---------------------------------------------------------------- physics
#slide[Real physics, computed][
  #grid(columns: (1fr, 1fr), gutter: 0.6cm,
  [
    #physics(title: "QED, Peskin and Schroeder")[
      - vacuum polarisation: Ward identity $C_g + q^2 C_(q q) = 0$ *exactly*
      - Peskin (7.91) reproduced to $10^(-15)$
      - running of $alpha$ from the leptons at one loop: $Delta alpha (M_Z^2) = 0.031419$
      - electron $g - 2$: #U, #F from the vertex graph, Dirac algebra in FORM, $F_2(0) = alpha \/ (2 pi)$ *exactly*
    ]
  ],
  [
    #physics(title: [$pi^0 -> gamma gamma$ from the quark triangle])[
      - trace with $gamma_5$ in FORM: $4 i m thin epsilon(k_1, k_2, mu, nu)$
      - triangle from its graph polynomials
      - width $7.79$ eV #h(0.3em) vs PrimEx-II $7.80 plus.minus 0.12$ eV
      - lifetime $8.35 times 10^(-17)$ s #h(0.3em) vs PDG $8.43 times 10^(-17)$ s
    ]
  ])
  #v(0.3em)
  Both are worked out in the notebooks `examples/peskin_examples.ipynb` and `examples/feynsage_walkthrough.ipynb`.
]

// ---------------------------------------------------------------- Chluba: Compton
#let chl = json("out/chluba_feynsage.json")
#let chir = json("out/chluba_ir.json")
#let chfc = json("out/chluba_ir_feyncalc.json")
#slide[From Chluba's thesis (2005): Compton scattering][
  #text(12pt, fill: ink2)[J. Chluba, _Spectral distortions of the cosmic microwave background_ (2005): Compton, double Compton and their infrared divergence. Every QED result there was recomputed twice, with feynsage + FORM and with FeynCalc.]
  #align(center)[#grid(columns: 2, gutter: 1.6cm, compton(channel: "s", scale: 0.75), compton(channel: "u", scale: 0.75))]
  #v(-0.3em)
  #result(title: "Both programs give, exactly")[
    #set text(13pt)
    $ sum |cal(M)|^2 = 8 e^4 [ -(u - m^2)/(s - m^2) - (s - m^2)/(u - m^2) + 4 m^2 (1/(s-m^2) + 1/(u - m^2)) + 4 m^4 (1/(s-m^2) + 1/(u - m^2))^2 ] $
    $ (dif sigma)/(dif Omega) = r_0^2/2 (omega'/omega)^2 (omega'/omega + omega/omega' - sin^2 theta) quad "(Klein–Nishina)", quad quad sigma_T = (8 pi)/3 r_0^2 quad "(Thomson limit)" $
  ]
  #text(11pt, fill: ink2)[FeynCalc: 0.2 s, symbolic. feynsage: exact rationals through FORM. Script: `examples/chluba/`.]
]

#slide[Double Compton: six diagrams][
  $e(P) + gamma(K_0) -> e(P') + gamma(K_1) + gamma(K_2)$: the photons attach to the electron line in all $3! = 6$ orders.
  #v(0.2em)
  #align(center)[#grid(columns: 3, gutter: 0.25cm,
    dc(("0", "1", "2"), scale: 0.6), dc(("0", "2", "1"), scale: 0.6), dc(("1", "0", "2"), scale: 0.6),
    dc(("2", "0", "1"), scale: 0.6), dc(("1", "2", "0"), scale: 0.6), dc(("2", "1", "0"), scale: 0.6))]
  #grid(columns: (1.25fr, 1fr), gutter: 0.6cm,
  [
    #result(title: "Against Mandl & Skyrme's X (thesis eq. D.1)")[
      #set text(13pt)
      At four random rational points both programs give *exactly* $sum_("spins, pol.") |cal(M)|^2 = 4 thin e^6 X$: the thesis's $|cal(M)|^2 = e^6 X$ is the average over the $2 times 2$ initial states.
    ]
  ],
  [
    #table(columns: 2, table.header([*36 traces, 12 $gamma$'s*], [*time per point*]),
      [feynsage + FORM], [0.5 s], [FeynCalc], [7.7 s])
  ])
]

#slide[The infrared divergence (thesis sec. 4.4.5)][
  #grid(columns: (0.85fr, 1.3fr), gutter: 0.5cm, align: horizon,
  [
    #align(center)[#soft-real(leg: "in", scale: 0.8) #h(0.2cm) #soft-real(leg: "out", scale: 0.8)]
    #align(center)[#text(10pt, fill: ink2)[a soft photon $K_2$ from an external electron]]
  ],
  [
    For $K_2 -> 0$ double Compton factorises into Compton times the *eikonal factor*:
    $ sum |cal(M)_"DC"|^2 -> e^2 S(K_2) sum |cal(M)_"C"|^2 $
    $ S = (2 P dot P')/((P dot K_2)(P' dot K_2)) - m^2/(P dot K_2)^2 - m^2/(P' dot K_2)^2 $
    #table(columns: 4, table.header([$K_2$ scaled by], [$10^(-4)$], [$10^(-6)$], [$10^(-8)$]),
      [ratio of the two sides], ..chl.soft.map(r => [#calc.round(r.at(1), digits: 10)]))
  ])
  $S prop 1\/omega_2^2$, so $dif sigma prop dif omega_2 \/ omega_2$: the divergence of the thesis. Integrated over angles, a soft photon is emitted with $dif N = alpha\/pi thin I(t) thin dif omega_2\/omega_2$, $I(t) = 2 P dot P' integral_0^1 (dif x)/(m^2 - x(1-x) t) - 2$.
  #physics(title: "Lightman's law (thesis eq. 4.24), derived")[
    For cold electrons and soft photons $I(t) = -2t\/3 + dots$; averaged over Thomson scattering, $-t = 2 omega_0^2 (1 - cos theta) -> 2 omega_0^2$, so
    $dif N \/ dif omega_2 = (4 alpha)/(3 pi) thin omega_0^2 \/ omega_2$ per Compton scattering: exactly the $4 alpha\/3 pi$ of the thesis.
  ]
]

#slide[Regularisation: the divergence cancels][
  #grid(columns: (0.75fr, 1.4fr), gutter: 0.5cm,
  [
    #align(center)[#soft-virtual(scale: 0.8)]
    #align(center)[#text(10pt, fill: ink2)[the virtual soft photon between the electron legs]]
    #v(0.2em)
    #align(center)[#vertex-loop(scale: 0.8)]
    #align(center)[#text(10pt, fill: ink2)[its pole is that of the QED vertex at $t = (P - P')^2$]]
  ],
  [
    The thesis stops the divergence with a lowest frequency $nu_(2,"min")$. In $D = 4 - 2 epsilon$ it cancels:
    $ "real" (omega_2 < Delta E): -(alpha I(t))/(2 pi epsilon), quad quad "virtual" (2 thin "Re" F_1): +(alpha I(t))/(2 pi epsilon) $
    The vertex, computed with the loop in $D$ dimensions, its IR pole from the soft triangle $C_0(m^2, t, m^2; 0, m, m)$ (UV poles do not depend on $t$, so differences isolate the IR part):
    #table(columns: 4,
      table.header([$t$], [$I(t) - I(-1\/2)$], [feynsage], [FeynCalc + Package-X]),
      ..range(3).map(i => {
        let r = chir.ir_lines.at(i); let f = chfc.lines.at(i)
        ([#r.at(0)], [#calc.round(r.at(2), digits: 12)], [#calc.round(r.at(1), digits: 12)], [#calc.round(f.at(1), digits: 12)])
      }).flatten())
    #text(11pt)[Normalisation check: $F_2 -> alpha\/2 pi$ as $t -> 0$ in both (1.99999967 at $t = -10^(-6)$, in units $alpha\/4 pi$).]
    #result(title: "So")[the cutoff $nu_(2,"min")$ of the thesis becomes $log(Delta E)$, the energy resolution: Bloch–Nordsieck, computed.]
  ])
]

// ---------------------------------------------------------------- comparison: features
#slide[Comparison with other packages][
  #set text(12pt)
  #table(columns: (25%, 11%, 11%, 11%, 11%, 10%, 10%, 11%),
    table.header([], [*feynsage*], [*Package-X*], [*LoopTools*], [*FeynCalc*], [*LiteRed*], [*FIRE*], [*Kira*]),
    [needs Mathematica], [no], [yes], [no], [yes], [yes], [yes#super[a]], [no],
    [one-loop tensor reduction], ok, ok, [numbers], ok, [], [], [],
    [$C_0$ as a formula], ok, ok, [], [via Package-X], [], [], [],
    [$D_0$ as a formula], ok, [IR cases], [numbers], [via Package-X], [], [], [],
    [IR poles of $C_0$, $D_0$], ok, ok, ok, [via Package-X], [], [], [],
    [multi-loop IBP], ok, [], [], [], ok, ok, ok,
    [finite fields], ok, [], [], [], [], ok, ok,
    [#U, #F from the graph], ok, [], [], ok, [], [], [],
    [Dirac algebra], [FORM], [`Spur`], [], ok, [], [], [],
  )
  #text(10pt, fill: ink2)[An empty cell: we did not find this feature in that program. #super[a]FIRE has a C++ core and a Mathematica front end.]
]

// ---------------------------------------------------------------- comparison: numbers
#slide[Head to head: the same job, measured][
  #grid(columns: (1.05fr, 1fr), gutter: 0.6cm,
  [
    #text(weight: "semibold")[Speed] (same laptop, symbols everywhere)
    #table(columns: (44%, 28%, 28%),
      table.header([*task*], [*feynsage*], [*Package-X*]),
      [$C_(001)$], [0.42 s], [0.95 s],
      [$C_(112)$], [0.39 s], [1.59 s],
      [$D_(00)$], [6.0 s], [0.08 s],
      [10 numerical $D_0$], [0.23 s], [0.25 s],
      [digits of a numerical $D_0$], [≈ 40], [16],
    )
    #text(10.5pt, fill: ink2)[Package-X is much faster for symbolic box tensors: honest numbers.]
  ],
  [
    #text(weight: "semibold")[Agreement]
    #table(columns: (50%, 50%),
      table.header([*against*], [*result*]),
      [Package-X, 190 points], [every digit of its values],
      [Package-X, IR, 118 points], [double precision],
      [LoopTools 2.16 (Denner), hard $D_0$ points], [$10^(-15)$],
      [Kira 3.1, LiteRed], [identical coefficients],
      [FeynCalc, #U and #F], [identical polynomials],
      [$D_0$ relabellings], [all 24 at 15 points; 6 random ones at 2000 random points],
    )
  ])
]

// ---------------------------------------------------------------- the bug
#slide[A point where Package-X and LoopTools fail][
  $C_0(-4, -3, -1; 0, 0, 1)$: every invariant is spacelike, so $C_0 = -integral dif^3 x thin delta(1 - sum x_i) \/ #F$ must be *real and negative*.
  #v(0.2em)
  #table(columns: (46%, 54%),
    table.header([*program*], [*value*]),
    [feynsage (closed form)], [#raw(read("out/bug_c0.txt").trim())],
    [Mathematica `NIntegrate` (no Package-X)], [`-0.585976809672364723`],
    [FeynCalc `FCFeynmanParametrize`, integrated], [`-0.58597680967236472`],
    [mpmath, 45 digits], [`-0.585976809672364722650390572218069267`],
    [Package-X 2.1.1, two labellings], text(fill: c-bad)[`+0.671253`, `-0.721045`],
    [LoopTools 2.16, version a (FF)], text(fill: c-bad)[`+0.671253`, `-0.721045`, `-1.843207`],
    [LoopTools 2.16, version b (Denner)], text(fill: c-bad)[four different complex numbers],
  )
  At nearby points Package-X and LoopTools' FF version agree with `NIntegrate` to 16 digits (the Denner version is still wrong there). Package-X's wrong values are LoopTools' FF values digit for digit.
]

// ---------------------------------------------------------------- tests
#slide[How we know it is right][
  #grid(columns: (1fr, 1fr), gutter: 0.6cm,
  [
    - 12 test files and a random stress test, run on a laptop (Sage 10.9) and on a CentOS 7 server (Sage 10.7)
    - every formula checked against another program *or* a direct integral
    - $D_0$ is symmetric under the 24 relabellings of its lines: all 24 agree at 15 points, 6 random ones at each of 2000 random points
    - exact expressions are evaluated at two precisions before they are returned
    - two independent code audits; every finding fixed and turned into a test
  ],
  [
    #note(title: "Limits, honestly")[
      - finite-field IBP: $D$ and at most one invariant
      - closed forms of $C_0$, $D_0$ need numerical kinematics; with symbols they stay symbols
      - vanishing Gram determinant: one group of identical lines, or two groups with equal masses
      - symbolic box tensors are slower than Package-X
    ]
  ])
]

// ---------------------------------------------------------------- end
#slide[Summary][
  #v(0.3em)
  #grid(columns: (1fr, 1fr), gutter: 0.8cm,
  [
    #result(title: "feynsage")[
      - graphs → #U, #F → families → IBP → masters
      - one loop: exact tensor reduction, $C_0$ and $D_0$ in closed form, IR poles, zero Gram determinants
      - QED of Chluba's thesis: Compton, double Compton, the IR divergence cancelled in $D$ dimensions
      - free, in SageMath, checked against Package-X, LoopTools, FeynCalc, LiteRed and Kira
    ]
  ],
  [
    ```bash
    git clone https://github.com/aburousan/feynsage
    cd feynsage && ./install.sh --test
    ```
    Start with `examples/feynsage_walkthrough.ipynb`.

    #v(0.4em)
    #text(fill: ink2)[Thanks to Prof. B. Ananthanarayan for the lecture and the `uf-new-short.nb` notebook that started the #U, #F part.]
  ])
]
