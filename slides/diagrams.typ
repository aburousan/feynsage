// Feynman diagrams for the slides, in the style of the lecture note (cetz; dark blue lines,
// photons as sine waves, fermions with an arrow in the middle).
#import "@preview/cetz:0.4.2": canvas, draw

#let line-col = rgb("#1f2d4d")
#let soft-col = rgb("#c8457a")
#let _vtx(pt) = draw.circle(pt, radius: 0.06, fill: line-col, stroke: none)
#let _lab(pt, body, size: 9pt, col: line-col) = draw.content(pt, text(size, fill: col, body))
#let _unit(a, b) = {
  let dx = b.at(0) - a.at(0); let dy = b.at(1) - a.at(1)
  let L = calc.sqrt(dx * dx + dy * dy)
  (dx / L, dy / L, L)
}
#let photon(a, b, amp: 0.09, waves: 5, col: line-col) = {
  let (ux, uy, L) = _unit(a, b)
  let n = 60
  let pts = range(n + 1).map(i => {
    let t = i / n
    let w = amp * calc.sin(t * waves * 2 * calc.pi)
    (a.at(0) + t * L * ux - w * uy, a.at(1) + t * L * uy + w * ux)
  })
  draw.line(..pts, stroke: (paint: col, thickness: 1.1pt))
}
#let fermion(a, b, col: line-col) = {
  draw.line(a, b, stroke: (paint: col, thickness: 1.25pt))
  let m1 = (a.at(0) * 0.56 + b.at(0) * 0.44, a.at(1) * 0.56 + b.at(1) * 0.44)
  let m2 = (a.at(0) * 0.44 + b.at(0) * 0.56, a.at(1) * 0.44 + b.at(1) * 0.56)
  draw.line(m1, m2, stroke: (paint: col, thickness: 1.25pt), mark: (end: "stealth", fill: col, scale: 0.8))
}
// photon arc between two points, bulging by h (for the virtual soft photon)
#let photon-arc(a, b, h: 0.8, waves: 7, amp: 0.08, col: line-col) = {
  let n = 90
  let pts = range(n + 1).map(i => {
    let t = i / n
    let x = a.at(0) + t * (b.at(0) - a.at(0))
    let y0 = a.at(1) + t * (b.at(1) - a.at(1)) + 4 * h * t * (1 - t)
    let w = amp * calc.sin(t * waves * 2 * calc.pi)
    (x, y0 + w)
  })
  draw.line(..pts, stroke: (paint: col, thickness: 1.1pt))
}

// ---------------------------------------------------------------- Compton
#let compton(channel: "s", scale: 1) = canvas(length: scale * 1cm, {
  let A = (-0.6, 0); let B = (0.6, 0)
  fermion((-1.7, -0.9), A); fermion(A, B); fermion(B, (1.7, -0.9))
  if channel == "s" {
    photon((-1.7, 0.9), A, waves: 4); photon(B, (1.7, 0.9), waves: 4)
    _lab((-1.95, 1.05), $K$); _lab((1.95, 1.05), $K'$)
  } else {
    photon((-1.7, 0.9), B, waves: 6); photon(A, (1.7, 0.9), waves: 6)
    _lab((-1.95, 1.05), $K$); _lab((1.95, 1.05), $K'$)
  }
  _vtx(A); _vtx(B)
  _lab((-1.95, -1.05), $P$); _lab((1.95, -1.05), $P'$)
})

// ---------------------------------------------------------------- double Compton: one order of the photons
// order: the photons attached from left to right along the electron line, e.g. ("0", "1", "2")
#let dc(order, scale: 1) = canvas(length: scale * 1cm, {
  let xs = (-0.9, 0, 0.9)
  fermion((-1.8, 0), (xs.at(0), 0)); fermion((xs.at(0), 0), (xs.at(1), 0))
  fermion((xs.at(1), 0), (xs.at(2), 0)); fermion((xs.at(2), 0), (1.8, 0))
  for (i, ph) in order.enumerate() {
    let v = (xs.at(i), 0)
    if ph == "0" {
      photon((xs.at(i) - 0.35, -1.0), v, waves: 3)
      _lab((xs.at(i) - 0.45, -1.25), $K_0$)
    } else {
      photon(v, (xs.at(i) + 0.3, 1.0), waves: 3)
      _lab((xs.at(i) + 0.42, 1.25), if ph == "1" { $K_1$ } else { $K_2$ })
    }
    _vtx(v)
  }
  _lab((-2.05, 0.2), $P$); _lab((2.05, 0.2), $P'$)
})

// ---------------------------------------------------------------- soft photon on an external leg
#let blob(c) = draw.circle(c, radius: 0.42, fill: luma(232), stroke: (paint: line-col, thickness: 1pt))
#let soft-real(leg: "out", scale: 1) = canvas(length: scale * 1cm, {
  let C = (0, 0)
  photon((-1.5, 0.9), (-0.35, 0.22), waves: 3); photon((0.35, 0.22), (1.5, 0.9), waves: 3)
  fermion((-1.6, -0.9), (-0.36, -0.22)); fermion((0.36, -0.22), (1.6, -0.9))
  blob(C)
  let v = if leg == "out" { (1.0, -0.6) } else { (-1.0, -0.6) }
  photon(v, (v.at(0) * 1.15, -1.7), waves: 3, col: soft-col)
  _vtx(v)
  _lab((v.at(0) * 1.15 + 0.35, -1.75), $K_2$, col: soft-col)
  _lab((-1.85, -1.0), $P$); _lab((1.85, -1.0), $P'$)
})
#let soft-virtual(scale: 1) = canvas(length: scale * 1cm, {
  let C = (0, 0)
  photon((-1.5, 0.9), (-0.35, 0.22), waves: 3); photon((0.35, 0.22), (1.5, 0.9), waves: 3)
  fermion((-1.6, -0.9), (-0.36, -0.22)); fermion((0.36, -0.22), (1.6, -0.9))
  blob(C)
  let a = (-1.1, -0.66); let b = (1.1, -0.66)
  photon-arc(a, b, h: -0.75, col: soft-col)
  _vtx(a); _vtx(b)
  _lab((-1.85, -1.0), $P$); _lab((1.85, -1.0), $P'$)
})
// the one-loop vertex that carries the same infrared pole
#let vertex-loop(scale: 1) = canvas(length: scale * 1cm, {
  let V = (0, 0.55); let A = (-0.75, -0.25); let B = (0.75, -0.25)
  photon((0, 1.6), V, waves: 3)
  fermion((-1.5, -1.05), A); fermion(A, V); fermion(V, B); fermion(B, (1.5, -1.05))
  photon(A, B, waves: 5, col: soft-col)
  for v in (V, A, B) { _vtx(v) }
  _lab((-1.75, -1.15), $P$); _lab((1.75, -1.15), $P'$); _lab((0.0, -0.6), $ell$, col: soft-col)
})
