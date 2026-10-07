r"""
Goncharov polylogarithms (GPLs) for solving canonical differential equations  df = eps dA f  order by order.

A function of two variables is kept as a sum  sum c * G(wx; x) G(wy; y)  with words wx, wy (tuples of
letters) and coefficients c that do not depend on x or y.  G(a1, ..., an; z) = Int_0^z dt/(t - a1) G(a2, ..., an; t)
with G(; z) = 1 and G(0, ..., 0; z) = log(z)^n/n!.  Integrating dt/(t - a) G(w; t) from 0 to z only puts the
letter a in front: G(a, w; z).  So with A = sum_k a_k log(letter_k), the order eps^n of f is obtained from
order eps^(n-1) by a matrix product and one new letter in every word.

    f1 = solve_canonical(fvec0, x_letters, y_letters, boundary)       # one order: straight path
                                                                        # (0,0) -> (x,0) -> (x,y)
    to_log(e)          weight <= 2 words with letters 0, 1, -1 in logarithms and dilogarithms
    at_one(e, table)   x = y = 1 with the values of G(w; 1) from a table (TABLE_AT_ONE by default)
"""
from sage.all import SR, log, polylog, pi, zeta, matrix


class GPL:
    """sum_k c_k G(wx_k; x) G(wy_k; y), stored as {(wx, wy): c}."""

    def __init__(self, terms=None):
        self.t = {}
        for k, c in (terms or {}).items():
            self._add(k, c)

    def _add(self, k, c):
        c = SR(c)
        if c.is_trivial_zero():
            return
        k = (tuple(SR(a) for a in k[0]), tuple(SR(a) for a in k[1]))
        key = (tuple(str(a) for a in k[0]), tuple(str(a) for a in k[1]))
        if key in self.t:
            old_k, old_c = self.t[key]
            nc = (old_c + c).simplify_full()
            if nc.is_trivial_zero():
                del self.t[key]
            else:
                self.t[key] = (old_k, nc)
        else:
            self.t[key] = (k, c)

    @staticmethod
    def const(c):
        return GPL({((), ()): c})

    @staticmethod
    def word(wx=(), wy=(), c=1):
        return GPL({(tuple(wx), tuple(wy)): c})

    def items(self):
        return [(k, c) for k, c in self.t.values()]

    def __add__(self, other):
        out = GPL()
        for k, c in self.items() + (other.items() if isinstance(other, GPL) else GPL.const(other).items()):
            out._add(k, c)
        return out

    __radd__ = __add__

    def __neg__(self):
        return GPL({k: -c for k, c in self.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rmul__(self, c):
        return GPL({k: SR(c) * v for k, v in self.items()})

    __mul__ = __rmul__

    def is_zero(self):
        return not self.t

    def prepend_x(self, a):
        """Int_0^x dt/(t - a) of this function at y = 0: words in y vanish there."""
        out = GPL()
        for (wx, wy), c in self.items():
            if wy:
                continue
            out._add(((SR(a),) + wx, ()), c)
        return out

    def prepend_y(self, a):
        out = GPL()
        for (wx, wy), c in self.items():
            out._add((wx, (SR(a),) + wy), c)
        return out

    def __repr__(self):
        if not self.t:
            return '0'
        parts = []
        for (wx, wy), c in self.items():
            g = ''
            if wx:
                g += 'G(%s; x)' % ', '.join(str(a) for a in wx)
            if wy:
                g += ('*' if g else '') + 'G(%s; y)' % ', '.join(str(a) for a in wy)
            parts.append('(%s)%s' % (c, ('*' + g) if g else ''))
        return ' + '.join(parts)

    def _latex_(self):
        from sage.all import latex
        if not self.t:
            return '0'
        parts = []
        for (wx, wy), c in self.items():
            g = ''
            if wx:
                g += r'G(%s;x)' % ','.join(latex(a) for a in wx)
            if wy:
                g += r'G(%s;y)' % ','.join(latex(a) for a in wy)
            cc = latex(c)
            parts.append((r'\left(%s\right) ' % cc if g else cc) + g if not (c - 1).is_trivial_zero() or not g else g)
        return ' + '.join(parts)

    def to_sr(self, rules):
        """Replace every G by rules(word, var) (a function returning an expression)."""
        tot = SR(0)
        x, y = SR.var('x'), SR.var('y')
        for (wx, wy), c in self.items():
            tot += c * rules(wx, x) * rules(wy, y)
        return tot


def solve_canonical(prev, x_part, y_part, boundary):
    r"""
    One order of  df = eps dA f.  prev: list of GPL (the previous order, one per integral).
    x_part: [(matrix, letter)] for the path (0, 0) -> (x, 0)  (the forms dlog(x - letter) at y = 0),
    y_part: [(matrix, letter)] for the path (x, 0) -> (x, y)  (the forms dlog(y - letter)),
    boundary: the constants at the starting point.  Returns the new list of GPL.
    """
    n = len(prev)
    out = [GPL.const(b) for b in boundary]
    y = SR.var('y')
    for M, a in x_part:
        a = SR(a).subs({y: 0})                       # the first leg runs along y = 0
        M = M.apply_map(lambda c: SR(c).subs({y: 0}))
        for i in range(n):
            acc = GPL()
            for j in range(n):
                if M[i, j] != 0:
                    acc = acc + M[i, j] * prev[j]
            out[i] = out[i] + acc.prepend_x(a)
    for M, a in y_part:
        for i in range(n):
            acc = GPL()
            for j in range(n):
                if M[i, j] != 0:
                    acc = acc + M[i, j] * prev[j]
            out[i] = out[i] + acc.prepend_y(a)
    return out


def to_log(word, z):
    """G(word; z) in logarithms and dilogarithms, for weight <= 2 and letters 0, 1, -1."""
    w = tuple(str(a) for a in word)
    table = {(): 1, ('0',): log(z), ('1',): log(1 - z), ('-1',): log(1 + z),
             ('0', '0'): log(z)**2/2, ('1', '0'): polylog(2, z) + log(z)*log(1 - z), ('0', '1'): -polylog(2, z),
             ('1', '1'): log(1 - z)**2/2, ('-1', '0'): polylog(2, -z) + log(z)*log(1 + z),
             ('0', '-1'): -polylog(2, -z), ('-1', '-1'): log(1 + z)**2/2}
    if w not in table:
        raise KeyError("no logarithm form for G(%s; z) here" % ', '.join(w))
    return table[w]


Lnybar = SR.var('Lnybar', latex_name=r'\ln\bar y')       # log(1 - y) at y = 1, where it diverges

TABLE_AT_ONE = {(): 1, ('0',): 0, ('1',): Lnybar, ('-1',): log(2), ('0', '0'): 0, ('1', '0'): pi**2/6,
                ('0', '1'): -pi**2/6, ('1', '1'): Lnybar**2/2, ('-1', '0'): -pi**2/12,
                ('0', '0', '0'): 0, ('1', '0', '0'): -zeta(3), ('-1', '0', '0'): 3*zeta(3)/4, ('0', '1', '0'): 2*zeta(3),
                ('0', '-1', '0'): -3*zeta(3)/2, ('0', '0', '1'): -zeta(3), ('0', '1', '1'): zeta(3),
                ('-1', '-1', '0'): -pi**2*log(2)/12 + zeta(3)/8, ('-1', '0', '1'): -pi**2*log(2)/6 + 5*zeta(3)/8,
                ('1', '1', '0'): Lnybar*pi**2/6 + zeta(3), ('1', '0', '1'): -Lnybar*pi**2/6 - 2*zeta(3),
                ('1', '1', '1'): Lnybar**3/6}


def at_one(e, table=None):
    """The value at x = y = 1 (letters that depend on x are evaluated at x = 1 first)."""
    table = table or TABLE_AT_ONE
    x = SR.var('x')
    tot = SR(0)
    for (wx, wy), c in e.items():
        kx = tuple(str(SR(a).subs({x: 1})) for a in wx)
        ky = tuple(str(SR(a).subs({x: 1})) for a in wy)
        tot += SR(c).subs({x: 1, SR.var('y'): 1}) * table[kx] * table[ky]
    return tot.simplify_full()
