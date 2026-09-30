r"""
Plots in the style of the lecture notes: the same pink and blue colours, Computer
Modern fonts, light grid, no top and right frame.

- set_theme()                  apply the style to matplotlib, and so to Sage's own plot()
- curves(funcs, (a, b), ...)   several real functions on one axis, from Sage expressions or callables
- mb_plane(left, right, c)     poles of a Mellin-Barnes integrand and the straight contour
- draw_graph(g, ...)           a FeynmanGraph as a Feynman diagram (also g.plot())

Every function returns a matplotlib Figure; save it with fig.savefig("name.svg").
"""
from itertools import permutations
from math import atan2, cos, pi, sin, sqrt

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

BLUE, PINK, VIOLET = "#2a78d6", "#d55181", "#4a3aa7"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
PALETTE = [BLUE, PINK, VIOLET, "#3a7fb5", "#b8336a"]

STYLE = {
    "font.family": "serif", "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm", "font.size": 10, "axes.formatter.use_mathtext": True,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 1.6,
    "axes.prop_cycle": mpl.cycler(color=PALETTE),
    "svg.fonttype": "path", "figure.dpi": 150, "axes.titlesize": 10,
    "legend.frameon": False,
}


def set_theme():
    """Use the style of the notes for every later matplotlib (and Sage) plot."""
    mpl.rcParams.update(STYLE)


def _fn(f):
    """A Sage expression in one variable, or any callable, as a float function."""
    if callable(f) and not hasattr(f, 'variables'):
        return lambda x: float(f(x))
    from sage.all import fast_callable, RDF
    v = f.variables()
    if len(v) != 1:
        raise ValueError("give functions of one variable")
    g = fast_callable(f, vars=v, domain=RDF)
    return lambda x: float(g(x))


def curves(funcs, xrange, labels=None, points=400, xlabel="$x$", ylabel=None, title=None,
           ylim=None, figsize=(4.2, 2.8)):
    r"""
    Plot several real functions on one axis.

        from feynsage.plotting import curves
        x = var('x')
        fig = curves([x/(x^2+0.1^2), 1/x], (-3, 3), labels=[r'$x/(x^2+\eta^2)$', '$1/x$'], ylim=(-6, 6))
        fig.savefig('pv.svg')
    """
    with mpl.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        xs = np.linspace(float(xrange[0]), float(xrange[1]), points)
        for k, f in enumerate(funcs):
            g = _fn(f)
            ys = []
            for x in xs:
                try:
                    ys.append(g(x))
                except (ZeroDivisionError, ValueError, OverflowError):
                    ys.append(np.nan)
            ys = np.array(ys)
            if ylim is not None:            # do not join the two sides of a pole
                ys[(ys < ylim[0] - (ylim[1] - ylim[0])) | (ys > ylim[1] + (ylim[1] - ylim[0]))] = np.nan
            ax.plot(xs, ys, color=PALETTE[k % len(PALETTE)], label=labels[k] if labels else None)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title, loc="left")
        if labels:
            ax.legend(fontsize=8)
        fig.tight_layout()
    return fig


def mb_plane(left, right, c, left_label=None, right_label=None, close=None, figsize=(5.0, 3.0)):
    r"""
    The complex z plane of a Mellin-Barnes integral: poles that must stay on the left
    (filled), poles that must stay on the right (open), and the contour Re z = c.
    close = "right" or "left" adds the arc used to close the contour.

        mb_plane([-1.4, -2.4, -3.4], [0, 1, 2, 3], -0.7,
                 r'$\Gamma(\lambda+z)$', r'$\Gamma(-z)$', close='right')
    """
    with mpl.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        ax.grid(False)
        for s in ["left", "bottom"]:
            ax.spines[s].set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
        xs = [float(v) for v in list(left) + list(right) + [c]]
        lo, hi = min(xs) - 1, max(xs) + 1.2
        H = 0.35 * (hi - lo)
        ax.annotate("", xy=(hi, 0), xytext=(lo, 0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=.8))
        ax.annotate("", xy=(0, H), xytext=(0, -H), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=.8))
        ax.text(hi, -.06 * H, r"Re $z$", color=INK2, ha="right", va="top")
        ax.text(.03 * H, .92 * H, r"Im $z$", color=INK2)
        ax.plot([float(v) for v in left], [0] * len(left), "o", color=BLUE, ms=7)
        ax.plot([float(v) for v in right], [0] * len(right), "o", mfc="white", mec=PINK, mew=1.6, ms=7)
        if left_label:
            ax.text(sum(map(float, left)) / len(left), .15 * H, "poles of " + left_label, color=BLUE, ha="center")
        if right_label:
            ax.text(sum(map(float, right)) / len(right), .15 * H, "poles of " + right_label, color=PINK, ha="center")
        c = float(c)
        ax.plot([c, c], [-.95 * H, .95 * H], color=VIOLET, lw=1.8, ls=(0, (5, 2)))
        for y in (-.6 * H, .4 * H):
            ax.annotate("", xy=(c, y + .15 * H), xytext=(c, y), arrowprops=dict(arrowstyle="-|>", color=VIOLET, lw=1.8))
        ax.text(c - .05, .85 * H, r"Re $z=c$", color=VIOLET, ha="right", fontsize=9)
        if close in ("right", "left"):
            sgn = 1 if close == "right" else -1
            width = (hi - c) if sgn > 0 else (c - lo)
            th = np.linspace(pi / 2, -pi / 2, 100)
            ax.plot(c + sgn * .85 * width * np.cos(th), .95 * H * np.sin(th), color=INK2, lw=1, ls=":")
        ax.set_xlim(lo - .2, hi + .2); ax.set_ylim(-H, 1.05 * H)
        ax.set_aspect("equal")
    return fig


# ---------------------------------------------------------------------- Feynman diagrams
def _crossings(order, edges):
    """Number of crossing chords when the vertices sit on a circle in this order."""
    pos = {v: i for i, v in enumerate(order)}
    chords = [tuple(sorted((pos[u], pos[v]))) for u, v in edges if u != v]
    n = 0
    for i in range(len(chords)):
        a, b = chords[i]
        for j in range(i + 1, len(chords)):
            c, d = chords[j]
            if a < c < b < d or c < a < d < b:
                n += 1
    return n


def _momentum_tex(q):
    parts = []
    for name, c in q.items():
        s = name if len(name) == 1 else "%s_{%s}" % (name[0], name[1:])
        if c == 1:
            parts.append("+" + s)
        elif c == -1:
            parts.append("-" + s)
        else:
            parts.append("%+d%s" % (c, s))
    t = "".join(parts)
    return "$" + (t[1:] if t.startswith("+") else t) + "$"


def draw_graph(g, labels=True, dots=None, momenta=True, figsize=(2.8, 2.8), R=1.0, removed=(), ax=None,
               title=None, names=None):
    r"""
    Draw a FeynmanGraph.  Vertices sit on a circle in the order with fewest crossing
    lines; external legs point outwards with their incoming momentum; massive lines are
    thick and blue, massless lines thin and black; line i is labelled x_i (as in U and F);
    dots = {i: n} puts n pink dots on line i (the propagator raised to n more powers);
    removed = lines (numbered from 1) drawn faint and dashed, as in the spanning-tree
    pictures of the notes; ax = draw into an existing matplotlib axis; names = the
    numbers to print on the lines (default 1, 2, ...), useful after contract().
    """
    dots = dots or {}
    removed = set(removed)
    V = list(g.vertices)
    edges = [(u, v) for u, v, _ in g.lines]
    # vertex order on the circle: fewest crossings, externals kept in their given order if possible
    ext_order = [v for v in g.ext]
    best, best_key = None, None
    cands = permutations(V[1:]) if len(V) <= 8 else [tuple(V[1:])]
    for rest in cands:
        order = (V[0],) + tuple(rest)
        pos = {v: i for i, v in enumerate(order)}
        key = (_crossings(order, edges), [pos[v] for v in ext_order] != sorted(pos[v] for v in ext_order))
        if best_key is None or key < best_key:
            best, best_key = order, key
    if ext_order:                                   # first external vertex on the left
        k0 = best.index(ext_order[0])
        best = best[k0:] + best[:k0]
    n = len(best)
    ang = {v: pi - 2 * pi * i / n if n > 2 else (pi if i == 0 else 0.0) for i, v in enumerate(best)}
    if n == 1:
        ang = {best[0]: 0.0}
    P = {v: (0.0, 0.0) if n == 1 else (R * cos(ang[v]), R * sin(ang[v])) for v in best}
    # a chain (bubbles in a row, sunsets): vertices on a horizontal line instead of a circle
    nbr = {v: set() for v in V}
    for u, w in edges:
        if u != w:
            nbr[u].add(w); nbr[w].add(u)
    if 2 < n and all(len(nbr[v]) <= 2 for v in V) and sum(len(x) for x in nbr.values()) == 2 * (n - 1):
        ends = [v for v in V if len(nbr[v]) == 1]
        start = ext_order[0] if ext_order and ext_order[0] in ends else ends[0]
        path, prev = [start], None
        while len(path) < n:
            nxt = [w for w in nbr[path[-1]] if w != prev][0]
            prev = path[-1]; path.append(nxt)
        best = tuple(path)
        P = {v: (-R + 2 * R * i / (n - 1), 0.0) for i, v in enumerate(path)}
        ang = {v: pi / 2 for v in path}
        ang[path[0]], ang[path[-1]] = pi, 0.0

    with mpl.rc_context(STYLE):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        ax.set_aspect("equal"); ax.axis("off")
        if title:
            ax.set_title(title, fontsize=9)
        groups = {}
        for i, (u, v, m2) in enumerate(g.lines):
            groups.setdefault(tuple(sorted((u, v), key=str)), []).append(i)
        for (u, v), idx in groups.items():
            k = len(idx)
            for s, i in enumerate(idx):
                m2 = g.lines[i][2]
                massive = (m2 != 0)
                col, lw = (BLUE, 2.6) if massive else (INK, 1.4)
                ls = "-"
                if (i + 1) in removed:
                    col, lw, ls = "#c9c7c2", 1.1, (0, (3, 2))
                if u == v:                                   # a tadpole loop on one vertex
                    (x0, y0), a = P[u], ang[u] + pi / 2 + s * 0.8
                    r = 0.45 * R
                    cx, cy = x0 + r * cos(a), y0 + r * sin(a)
                    ax.add_patch(Circle((cx, cy), r, fill=False, color=col, lw=lw, ls=ls))
                    lx, ly, mid = cx + 1.45 * r * cos(a), cy + 1.45 * r * sin(a), (cx + r * cos(a), cy + r * sin(a))
                    tang = a + pi / 2
                else:
                    (x1, y1), (x2, y2) = P[u], P[v]
                    L = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
                    nx, ny = -(y2 - y1) / L, (x2 - x1) / L
                    bend = (s - (k - 1) / 2) * 0.55
                    mx, my = (x1 + x2) / 2 + bend * L * nx / 2, (y1 + y2) / 2 + bend * L * ny / 2
                    # quadratic Bezier through (mx, my): control point 2*m - (p1+p2)/2
                    cx, cy = 2 * mx - (x1 + x2) / 2, 2 * my - (y1 + y2) / 2
                    t = np.linspace(0, 1, 60)
                    bx = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cx + t ** 2 * x2
                    by = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cy + t ** 2 * y2
                    ax.plot(bx, by, color=col, lw=lw, ls=ls, solid_capstyle="round")
                    side = 1 if bend >= 0 else -1
                    if k == 1:                               # label on the outer side of a single line
                        side = 1 if (mx * nx + my * ny) >= 0 else -1
                    lx, ly = mx + side * 0.22 * R * nx, my + side * 0.22 * R * ny
                    mid = (mx, my)
                    tang = atan2(y2 - y1, x2 - x1)
                if labels:
                    ax.text(lx, ly, "$x_{%d}$" % (names[i] if names else i + 1), ha="center", va="center", fontsize=9,
                            color="#b5b3ae" if (i + 1) in removed else (BLUE if massive else INK2))
                nd = dots.get(names[i] if names else i + 1, 0)
                for j in range(nd):
                    off = (j - (nd - 1) / 2) * 0.16 * R
                    ax.plot([mid[0] + off * cos(tang)], [mid[1] + off * sin(tang)], "o", color=PINK, ms=6, zorder=5)
        # external legs
        for vtx, q in g.ext.items():
            x0, y0 = P[vtx]
            a = ang[vtx] if n > 1 else pi
            x1, y1 = x0 + 0.62 * R * cos(a), y0 + 0.62 * R * sin(a)
            ax.plot([x0, x1], [y0, y1], color=INK, lw=1.4)
            ax.add_patch(FancyArrowPatch((x0 + .45 * R * cos(a), y0 + .45 * R * sin(a)),
                                         (x0 + .25 * R * cos(a), y0 + .25 * R * sin(a)),
                                         arrowstyle="-|>", mutation_scale=10, color=INK, lw=1))
            if momenta:
                ca = cos(a)
                ax.text(x0 + .7 * R * cos(a), y0 + .7 * R * sin(a) + (.12 * R * (1 if sin(a) >= 0 else -1) if abs(ca) <= .3 else 0),
                        _momentum_tex(q), ha="right" if ca < -.3 else ("left" if ca > .3 else "center"),
                        va="center", fontsize=10)
        for v in best:
            ax.plot(*P[v], "o", color=INK, ms=4, zorder=6)
        ax.set_xlim(-1.9 * R, 1.9 * R); ax.set_ylim(-1.9 * R, 1.9 * R)
    return fig


def draw_panels(g, removed_sets, titles=None, ncols=5, size=1.9, **kw):
    """Several copies of the diagram side by side, each with its own removed lines
    (for example all spanning trees, or all 2-forests)."""
    n = len(removed_sets)
    nrows = (n + ncols - 1) // ncols
    with mpl.rc_context(STYLE):
        fig, axes = plt.subplots(nrows, min(n, ncols), figsize=(size * min(n, ncols), size * nrows), squeeze=False)
        for k, ax in enumerate(axes.flat):
            if k < n:
                draw_graph(g, removed=[i + 1 for i in removed_sets[k]], ax=ax,
                           title=titles[k] if titles else None, **kw)
            else:
                ax.axis("off")
        fig.tight_layout()
    return fig


def draw_sectors(g, sectors, titles=None, ncols=5, size=1.9, **kw):
    """Sectors of a diagram: sector s (a 0/1 tuple, one entry per line) keeps the lines
    with 1 and shrinks the lines with 0 to points (see FeynmanGraph.contract)."""
    n = len(sectors)
    nrows = (n + ncols - 1) // ncols
    with mpl.rc_context(STYLE):
        fig, axes = plt.subplots(nrows, min(n, ncols), figsize=(size * min(n, ncols), size * nrows), squeeze=False)
        for k, ax in enumerate(axes.flat):
            if k < n:
                h, kept = g.contract([i + 1 for i, e in enumerate(sectors[k]) if not e])
                draw_graph(h, ax=ax, names=kept, title=titles[k] if titles else str(tuple(sectors[k])), **kw)
            else:
                ax.axis("off")
        fig.tight_layout()
    return fig
