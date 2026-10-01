# Benchmark families: the same integral family written for feynsage and for Kira 3.1.
# Each case: loops, propagators (momentum, mass^2) in Kira syntax, the external momenta, the scalar
# products, the symbol Kira replaces by one (feynsage gets the number 1 directly), the top sector
# (lines that are propagators; the rest are numerator lines) and the sizes to run.
CASES = {
    "kite": dict(
        loops=["l1", "l2"],
        props=[("l1", "0"), ("l1+q", "0"), ("l1+l2", "0"), ("l1+l2+q", "0"), ("l2", "0")],
        ext=["q"],
        sp={("q", "q"): "qq"},
        invariants=[("qq", 2)], one="qq",
        top=5,
        sizes=[(2, 1), (4, 2), (6, 2), (8, 3), (10, 4)],          # (extra dots r, numerator powers s)
    ),
    "vertex": dict(
        loops=["l1", "l2"],
        props=[("l1", "0"), ("l2", "0"), ("-l2+p1+p2", "0"), ("-l1+p1+p2", "0"), ("l1-p1", "0"),
               ("-l1+l2", "0"), ("l2+p1", "0")],
        ext=["p1", "p2"],
        sp={("p1", "p1"): "0", ("p2", "p2"): "0", ("p1", "p2"): "qq/2"},
        invariants=[("qq", 2)], one="qq",
        top=6,
        sizes=[(1, 1), (2, 2), (3, 2), (4, 3), (5, 3)],
    ),
    "sunset": dict(
        loops=["k", "l"],
        props=[("k", "msq"), ("k-l", "msq"), ("l+p", "msq"), ("k+p", "0"), ("l", "0")],
        ext=["p"],
        sp={("p", "p"): "pp"},
        invariants=[("pp", 2), ("msq", 2)], one="msq",
        top=3,
        sizes=[(2, 1), (4, 2), (6, 2), (8, 3), (10, 3)],
    ),
    "doublebox": dict(
        loops=["k1", "k2"],
        props=[("k1", "0"), ("k1+p1", "0"), ("k1+p1+p2", "0"), ("k2+p1+p2", "0"), ("k2+p1+p2+p3", "0"),
               ("k2", "0"), ("k1-k2", "0"), ("k1+p1+p2+p3", "0"), ("k2+p1", "0")],
        ext=["p1", "p2", "p3"],
        sp={("p1", "p1"): "0", ("p2", "p2"): "0", ("p3", "p3"): "0",
            ("p1", "p2"): "s/2", ("p2", "p3"): "t/2", ("p1", "p3"): "-s/2-t/2"},
        invariants=[("s", 2), ("t", 2)], one="t",
        top=7,
        sizes=[(0, 1), (1, 1), (1, 2), (2, 2), (3, 2)],
    ),
}


def targets(case, r, s):
    """All integrals of the top sector with exactly r extra dots and s numerator powers."""
    from itertools import product
    c = CASES[case]
    t, top = len(c["props"]), c["top"]

    def comps(n, k):
        if k == 0:
            if n == 0:
                yield ()
            return
        for first in range(n + 1):
            for rest in comps(n - first, k - 1):
                yield (first,) + rest
    out = []
    for dots in comps(r, top):
        for nums in (comps(s, t - top) if t > top else [()]):     # no numerator lines: s only sets the seeds
            out.append(tuple(1 + d for d in dots) + tuple(-n for n in nums))
    return out
