// Slide style, with the colours, fonts and dashed boxes of the lecture note
// "Feynman integrals at one loop".
#import "@preview/mitex:0.2.7": mitex

#let ink = rgb("#0b0b0b")
#let ink2 = rgb("#52514e")
#let c-blue = rgb("#2a64c4")
#let c-pink = rgb("#c8457a")
#let c-green = rgb("#3a7fb5")
#let c-violet = rgb("#4a3aa7")
#let c-grey = rgb("#6a6f86")
#let c-good = rgb("#2f8a57")
#let c-bad = rgb("#c0392b")

#let setup(body) = {
  set document(title: "feynsage", author: "Kazi Abu Rousan")
  set page(width: 25.4cm, height: 14.29cm, margin: (x: 1.3cm, top: 1.25cm, bottom: 0.9cm), fill: rgb("#fffdfd"),
    footer: context {
      set text(7.5pt, fill: ink2)
      [feynsage: Feynman integrals in SageMath #h(1fr) #counter(page).display()]
    })
  set text(font: "New Computer Modern", size: 15pt, fill: ink, lang: "en")
  show math.equation: set text(font: "New Computer Modern Math")
  set par(justify: false, leading: 0.55em, spacing: 0.75em)
  show raw: set text(font: "Menlo", size: 10.5pt)
  show raw.where(block: true): it => block(width: 100%, inset: (x: 8pt, y: 7pt), radius: 4pt,
    fill: c-violet.lighten(95%), stroke: (paint: c-violet.lighten(50%), thickness: 0.6pt, dash: "dashed"), it)
  set table(inset: 5pt, align: left + horizon,
    stroke: (x, y) => if y == 0 { (bottom: 0.8pt + ink2) } else { (bottom: 0.3pt + ink2.lighten(60%)) })
  show table: set text(size: 11.5pt)
  set list(indent: 0.4em, spacing: 0.55em, marker: text(fill: c-pink)[•])
  body
}

#let slide(title, body) = {
  pagebreak(weak: true)
  block(below: 0.55em)[
    #text(22pt, weight: "semibold", fill: c-blue.darken(25%))[#title]
    #v(-0.55em)
    #line(length: 100%, stroke: (paint: c-pink.lighten(20%), thickness: 0.9pt, dash: "dashed"))
  ]
  body
}

#let dbox(title, colour, body) = block(width: 100%, inset: (x: 9pt, y: 7pt), radius: 5pt,
  stroke: (paint: colour, thickness: 0.9pt, dash: "dashed"), fill: colour.lighten(96%))[
  #if title != none [#text(weight: "semibold", fill: colour.darken(15%), size: 10.5pt)[#upper(title)] \ ]
  #body
]
#let result(title: "Result", body) = dbox(title, c-blue, body)
#let physics(title: "Physical picture", body) = dbox(title, c-green, body)
#let derivation(title: "Derivation", body) = dbox(title, c-pink, body)
#let note(title: "Remark", body) = dbox(title, c-grey, body)
#let out(path, size: 10.5pt) = block(width: 100%, inset: (x: 8pt, y: 6pt), radius: 4pt,
  fill: c-blue.lighten(95%), stroke: (paint: c-blue.lighten(45%), thickness: 0.6pt))[
  #set text(font: "Menlo", size: size)
  #read(path).trim()
]
#let ok = text(fill: c-good, weight: "bold")[✓]
#let bad = text(fill: c-bad, weight: "bold")[✗]
#let U = $cal(U)$
#let F = $cal(F)$

// a Sage result, typeset from the LaTeX that Sage writes (latex(expr))
#let mt(path, size: 13pt) = block(width: 100%, inset: (x: 8pt, y: 6pt), radius: 4pt,
  fill: c-blue.lighten(95%), stroke: (paint: c-blue.lighten(45%), thickness: 0.6pt))[
  #set text(size: size)
  #mitex(read(path))
]
