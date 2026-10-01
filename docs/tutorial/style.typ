// Page style and boxes of the lecture note "Feynman integrals at one loop", with notebook cells.

#let ink = rgb("#0b0b0b")
#let ink2 = rgb("#52514e")
#let c-blue = rgb("#2a64c4")     // headings, results
#let c-pink = rgb("#c8457a")     // accents, derivations
#let c-orange = rgb("#b8336a")   // board corrections (deep rose)
#let c-green = rgb("#3a7fb5")    // physical pictures (soft steel blue)
#let c-violet = rgb("#4a3aa7")   // Mathematica
#let c-grey = rgb("#6a6f86")     // remarks

#let setup(body) = {
  set document(title: "feynsage: a step-by-step tutorial", author: "Kazi Abu Rousan")
  set page(
    paper: "a4",
    margin: (x: 1.55cm, top: 2.0cm, bottom: 1.9cm),
    numbering: "1",
    header: context {
      if counter(page).get().first() > 1 [
        #set text(8.5pt, fill: ink2)
        feynsage: a step-by-step tutorial #h(1fr) Feynman integrals in SageMath
        #v(-4pt)
        #line(length: 100%, stroke: (paint: c-pink.lighten(40%), thickness: 0.5pt, dash: "dashed"))
      ]
    },
  )
  set text(font: "New Computer Modern", size: 10.5pt, fill: ink, lang: "en")
  set page(fill: rgb("#fffdfd"))
  show math.equation: set text(font: "New Computer Modern Math")
  set par(justify: true, leading: 0.62em, spacing: 0.95em)
  set heading(numbering: "1.1")
  show heading.where(level: 1): it => {
    v(1.2em)
    block(below: 1em)[
      #set text(17pt, weight: "semibold", fill: c-blue.darken(25%))
      #if it.numbering != none [#counter(heading).display() #h(0.4em)]
      #it.body
      #v(-0.5em)
      #line(length: 100%, stroke: (paint: c-pink.lighten(20%), thickness: 0.9pt, dash: "dashed"))
    ]
  }
  show heading.where(level: 2): set text(12.5pt, weight: "semibold", fill: c-pink.darken(15%))
  show heading.where(level: 3): set text(11pt, weight: "semibold", fill: c-blue.darken(20%))
  show heading: set block(sticky: true)
  show raw: set text(font: "Menlo", size: 8.2pt)
  show link: set text(fill: c-blue)
  show table: set par(justify: false)
  show table: set text(size: 9.5pt)
  set table(align: left + horizon)
  set list(indent: 0.6em)
  set enum(indent: 0.6em)
  body
}

// generic dashed box
#let dbox(title, colour, body, breakable: true) = block(
  width: 100%,
  inset: (x: 10pt, y: 9pt),
  radius: 5pt,
  breakable: breakable,
  stroke: (paint: colour, thickness: 0.9pt, dash: "dashed"),
  fill: colour.lighten(97%),
  above: 1.1em,
  below: 1.1em,
)[
  #block(sticky: true, below: 0.5em)[#text(weight: "semibold", fill: colour.darken(15%), size: 9.5pt)[#upper(title)]]
  #body
]

#let result(title: "Result", body) = dbox(title, c-blue, body)
#let physics(title: "Physical picture", body) = dbox(title, c-green, body)
#let board(title: "Correction to the board", body) = dbox(title, c-orange, body)
#let derivation(title: "Derivation, step by step", body) = dbox(title, c-pink, body)
#let note(title: "Remark", body) = dbox(title, c-grey, body)
// a feynsage (or Sage) function, explained where it is first used
#let func(sig, body, what: "feynsage") = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 5pt, breakable: true,
  stroke: (paint: c-violet, thickness: 0.9pt, dash: "dashed"), fill: c-violet.lighten(96%), above: 1em, below: 1em)[
  #block(sticky: true, below: 0.55em)[#text(weight: "semibold", fill: c-violet.darken(10%), size: 9pt)[#upper(what) FUNCTION] #h(0.5em) #text(size: 9pt)[#raw(sig, lang: "python")]]
  #set par(justify: false)
  #body
]


// ------------------------------------------------------------------ notebook cells
// Every cell was run by run_cells.sage; out/cells.json holds its code, what it printed, the value of
// its last line (LaTeX or text) and its figure.
#import "@preview/mitex:0.2.7": mitex
#let _cells = json("out/cells.json")
#let _prompt(body, col) = text(7.5pt, fill: col, font: "Menlo")[#body]
#let _codebox(body) = block(width: 100%, inset: (x: 7pt, y: 5pt), radius: 3pt, fill: rgb("#f5f6fa"),
  stroke: 0.4pt + c-violet.lighten(70%))[#set text(font: "Menlo", size: 7.6pt); #set par(justify: false); #body]
#let _outbox(body) = block(width: 100%, inset: (x: 7pt, y: 4pt), radius: 3pt, fill: white,
  stroke: (left: 1.6pt + c-pink.lighten(45%)))[#set par(justify: false); #body]
#let cell(id, fig: 70%, math-size: 9.5pt) = {
  let c = _cells.at(id)
  let rows = (pad(top: 5pt, _prompt([In [#c.n]:], rgb("#3c5fa8"))), _codebox(raw(c.code, lang: "python", block: true)))
  if c.stdout != "" {
    rows += ([], _outbox[#set text(font: "Menlo", size: 7.4pt); #raw(c.stdout, block: true)])
  }
  if c.latex != none {
    rows += (pad(top: 4pt, _prompt([Out[#c.n]:], c-pink.darken(10%))), _outbox[#set text(size: math-size); #mitex(c.latex)])
  }
  if c.text != none {
    rows += (pad(top: 4pt, _prompt([Out[#c.n]:], c-pink.darken(10%))), _outbox[#set text(font: "Menlo", size: 7.4pt); #raw(c.text, block: true)])
  }
  if c.fig != none {
    rows += (pad(top: 4pt, _prompt([Out[#c.n]:], c-pink.darken(10%))), _outbox(align(center, image(c.fig, width: fig))))
  }
  block(width: 100%, above: 0.9em, below: 1em, breakable: true,
    grid(columns: (36pt, 1fr), column-gutter: 4pt, row-gutter: 4pt, align: (right + top, left + top), ..rows))
}
