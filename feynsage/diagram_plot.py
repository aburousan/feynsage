r"""
Drawing of the diagrams made by insert_fields (like FeynArts' Paint):

    diags.draw()              # a grid of all diagrams, a matplotlib Figure

Incoming legs on the left, outgoing legs on the right, internal vertices placed by a small spring
layout.  Fermions are solid lines with an arrow along the fermion-number flow, photons, W and Z
bosons wavy lines, gluons curly lines, scalars dashed and ghosts dotted.  Each line carries the
LaTeX name of its field.
"""
import math
import numpy as np


def _layout(top):
    """Positions of the end points: legs fixed on two sides, vertices by spring relaxation."""
    n_in, n_out = top.n_in, top.n_out
    pos = {}
    for i in range(n_in):
        pos[('x', i + 1)] = np.array([0.0, 1.0 - (i + 0.5) / n_in])
    for j in range(n_out):
        pos[('x', n_in + j + 1)] = np.array([1.0, 1.0 - (j + 0.5) / n_out])
    V = top.V
    rng = np.random.RandomState(7)
    for v in range(V):
        # start each vertex at the mean of its external legs (or the centre)
        legs = [pos[a] for (a, b) in top.lines if b == ('v', v) and a[0] == 'x'] + \
               [pos[b] for (a, b) in top.lines if a == ('v', v) and b[0] == 'x']
        start = np.mean(legs, axis=0) * 0.5 + np.array([0.25, 0.25]) if legs else np.array([0.5, 0.5])
        pos[('v', v)] = start + 0.02 * rng.randn(2)
    edges = [(a, b) for a, b in top.lines if a != b]
    for it in range(400):
        force = {('v', v): np.zeros(2) for v in range(V)}
        for a, b in edges:                                  # springs of length L0
            d = pos[b] - pos[a]
            r = np.linalg.norm(d) + 1e-9
            L0 = 0.28
            f = (r - L0) * d / r
            if a[0] == 'v':
                force[a] += f
            if b[0] == 'v':
                force[b] -= f
        for v in range(V):                                  # repulsion between vertices
            for w in range(V):
                if v == w:
                    continue
                d = pos[('v', v)] - pos[('v', w)]
                r = np.linalg.norm(d) + 1e-3
                force[('v', v)] += 0.004 * d / r ** 3
        for v in range(V):
            pos[('v', v)] = np.clip(pos[('v', v)] + 0.1 * force[('v', v)], 0.12, 0.88)
    return pos


def _wavy_unused(ax, P, Q, amp=0.018, waves=None, curly=False, **kw):
    d = Q - P
    L = np.linalg.norm(d)
    if L < 1e-9:
        return
    n = np.array([-d[1], d[0]]) / L
    waves = waves or max(3, int(L / 0.045))
    t = np.linspace(0, 1, 40 * waves)
    if curly:
        # a curly (gluon) line: loops made by a cycloid along the line
        ph = 2 * math.pi * waves * t
        x = P[0] + d[0] * t + (amp * 0.9) * (np.cos(ph) - 1) * d[0] / L * 0.6 + amp * np.sin(ph) * n[0]
        y = P[1] + d[1] * t + (amp * 0.9) * (np.cos(ph) - 1) * d[1] / L * 0.6 + amp * np.sin(ph) * n[1]
    else:
        off = amp * np.sin(2 * math.pi * waves * t)
        x = P[0] + d[0] * t + off * n[0]
        y = P[1] + d[1] * t + off * n[1]
    ax.plot(x, y, **kw)


def _path(P, Q, bend, npts=200):
    """Points on the straight line (bend = 0) or a circular-looking arc from P to Q."""
    t = np.linspace(0, 1, npts)
    d = Q - P
    L = np.linalg.norm(d) + 1e-9
    n = np.array([-d[1], d[0]]) / L
    C = (P + Q) / 2 + bend * L * n
    return np.outer((1 - t) ** 2, P) + np.outer(2 * (1 - t) * t, C) + np.outer(t ** 2, Q)


def _decorate(pts, kind, curly, amp=0.016):
    """Wavy or curly offset along a polyline."""
    seg = np.diff(pts, axis=0)
    length = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    tang = np.gradient(pts, axis=0)
    tang /= (np.linalg.norm(tang, axis=1)[:, None] + 1e-12)
    nrm = np.stack([-tang[:, 1], tang[:, 0]], axis=1)
    waves = max(3, int(length[-1] / 0.045))
    ph = 2 * math.pi * waves * length / (length[-1] + 1e-12)
    if curly:
        return pts + amp * np.sin(ph)[:, None] * nrm + 0.6 * amp * (np.cos(ph) - 1)[:, None] * tang
    return pts + amp * np.sin(ph)[:, None] * nrm


def _line(ax, P, Q, particle, label, arrow_forward, bend=0.0):
    kind = particle.kind
    col = 'black'
    pts = _path(P, Q, bend)
    if kind == 'V':
        w = _decorate(pts, kind, particle.colour == '8')
        ax.plot(w[:, 0], w[:, 1], color=col, lw=1.0)
    else:
        ls = '-' if kind == 'F' else ((0, (4, 3)) if kind == 'S' else (0, (1, 2)))
        ax.plot(pts[:, 0], pts[:, 1], color=col, lw=1.15, ls=ls)
    mid = pts[len(pts) // 2]
    if kind in ('F', 'U'):
        d = pts[len(pts) // 2 + 3] - pts[len(pts) // 2 - 3]
        d = d / (np.linalg.norm(d) + 1e-9) * (1 if arrow_forward else -1)
        ax.annotate('', xy=mid + 0.03 * d, xytext=mid - 0.03 * d,
                    arrowprops=dict(arrowstyle='-|>', color=col, lw=1.0, mutation_scale=9))
    if label:
        d = Q - P
        L = np.linalg.norm(d) + 1e-9
        n = np.array([-d[1], d[0]]) / L
        m = mid + (0.06 if bend >= 0 else -0.06) * n
        ax.text(m[0], m[1], '$%s$' % label, ha='center', va='center', fontsize=8, color='#2a64c4')


def _self_loop(ax, P, particle, label):
    th = np.linspace(-math.pi / 2, 3 * math.pi / 2, 160)
    r = 0.08
    pts = np.stack([P[0] + r * np.cos(th), P[1] + r + r * np.sin(th)], axis=1)
    if particle.kind == 'V':
        w = _decorate(pts, 'V', particle.colour == '8', amp=0.012)
        ax.plot(w[:, 0], w[:, 1], color='black', lw=1.0)
    else:
        ls = '-' if particle.kind == 'F' else ((0, (4, 3)) if particle.kind == 'S' else (0, (1, 2)))
        ax.plot(pts[:, 0], pts[:, 1], color='black', lw=1.1, ls=ls)
    if label:
        ax.text(P[0], P[1] + 2 * r + 0.04, '$%s$' % label, ha='center', fontsize=8, color='#2a64c4')


def draw_diagram(diag, ax):
    top, model = diag.top, diag.model
    pos = _layout(top)
    # lines joining the same two end points are bent apart
    groups = {}
    for k, (a, b) in enumerate(top.lines):
        if a != b:
            groups.setdefault(frozenset((a, b)), []).append(k)
    bends = {}
    for key, ks in groups.items():
        n = len(ks)
        for j, k in enumerate(ks):
            bends[k] = 0.0 if n == 1 else (j - (n - 1) / 2) * 0.8
    for k, (name, s, t) in enumerate(diag.fields):
        p = model.particles[name]
        a, b = top.lines[k]
        if a == b:
            _self_loop(ax, pos[a], p, p.latex)
            continue
        P, Q = pos[s], pos[t]
        bend = bends.get(k, 0.0)
        if (s, t) != tuple(sorted((a, b))):
            bend = -bend                              # bends are defined for one fixed order of the ends
        # the arrow follows the particle (not the antiparticle)
        forward = not name.endswith('~')
        lab = p.latex if not name.endswith('~') else model.particles[p.anti].latex
        if p.kind in ('F', 'U') and name.endswith('~'):
            lab = model.particles[p.anti].latex
        elif p.kind not in ('F', 'U'):
            lab = p.latex
        _line(ax, P, Q, p, lab, forward, bend)
    for v in range(top.V):
        ax.plot(*pos[('v', v)], 'o', color='black', ms=2.5)
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, 1.12)
    ax.set_aspect('equal')
    ax.axis('off')


def draw_diagrams(diags, ncols=4, size=2.2):
    import matplotlib.pyplot as plt
    n = max(1, len(diags))
    ncols = min(ncols, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(size * ncols, size * nrows), squeeze=False)
    for i, ax in enumerate(axes.flat):
        if i < len(diags):
            draw_diagram(diags[i], ax)
            ax.set_title('(%d)' % (i + 1), fontsize=8)
        else:
            ax.axis('off')
    fig.tight_layout()
    return fig
