r"""
FORM output -> Sage (internal).

FORM prints a fully expanded sum of monomials.  Each monomial is split into its factors, every
distinct factor (p.q, p(mu), d_(mu,nu), e_(...), a symbol, a power) is turned into a Sage object once
and cached, the monomial is the product of its factors and the sum is built with one n-ary addition.
This is about 15 times faster than eval() of the printed text (4383 terms: 23 ms against 320 ms).
Dirac matrices of open lines (g_, gi_, g5_, g6_, g7_) are returned as atom lists per line.
"""
import re
from sage.all import SR, QQ, I, Integer
from sage.misc.misc_c import prod
from . import tensors as T


class ParseError(ValueError):
    pass


def expression_text(output, name):
    """The printed body of expression `name` in FORM's output (whitespace and continuations removed)."""
    m = re.search(r'(?m)^\s*%s\s*=' % re.escape(name), output)
    if not m:
        raise ParseError("expression %s not found in FORM's output" % name)
    rest = output[m.end():]
    end = rest.find(';')
    body = rest[:end] if end >= 0 else rest
    body = body.replace('\\\n', '')
    return re.sub(r'\s+', '', body)


def _split_top(s, seps):
    """Split s at characters in seps at bracket depth 0; for '+-' a sign after '^' does not split."""
    out, depth, start = [], 0, 0
    for i, ch in enumerate(s):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif depth == 0 and ch in seps and i > start:
            if seps == '+-':
                if s[i - 1] in '^*/(,':
                    continue
                out.append(s[start:i])
                start = i                  # the sign stays with the next term
            else:
                out.append(s[start:i])
                start = i + 1
    out.append(s[start:])
    return out


_TERM = re.compile(r'(?<=[^\^*/(,])(?=[+-])')     # before a sign that starts a term
_ARG_OPS = re.compile(r'\([^()]*[+*]|\([^()]*(?<!\^)-')


def _args(s):
    """'name(a,b,c)' -> ('name', ['a', 'b', 'c']) for a flat argument list."""
    j = s.index('(')
    return s[:j], s[j + 1:-1].split(',')


class Reader:
    def __init__(self, program):
        self.p = program
        self.cache = {}
        self.dummies = {}

    # -- names back to Sage
    def _obj(self, n):
        hit = self.p.sage_of.get(n)
        if hit is not None:
            return hit[1]
        if re.fullmatch(r'N\d+_\?', n):                     # a dummy index made by FORM
            if n not in self.dummies:
                self.dummies[n] = T.lorentz_indices('dummy%d' % (len(self.dummies) + 1))
            return self.dummies[n]
        raise ParseError("FORM returned the unknown name %s" % n)

    def _slot_obj(self, n):
        return self._obj(n)

    def atom(self, s):
        """One factor of a monomial -> ('s', Sage scalar) or ('g', line, [Dirac atoms])."""
        hit = self.cache.get(s)
        if hit is not None:
            return hit
        r = self._atom(s)
        self.cache[s] = r
        return r

    def _atom(self, s):
        if s[0].isdigit():
            return ('s', SR(QQ(s)))
        if s == 'i_':
            return ('s', I)
        k = s.rfind('^')
        if k > 0 and s.count('(', k) == s.count(')', k) and s[k + 1:].lstrip('(-').rstrip(')').isdigit() and _depth_ok(s[:k]):
            base = self.atom(s[:k])
            if base[0] != 's':
                raise ParseError("power of a Dirac matrix in FORM output: %s" % s)
            n = int(s[k + 1:].strip('()'))
            return ('s', base[1] ** n)
        if '(' in s:
            name, args = _args(s)
            if name == 'd_':
                a, b = [self._slot_obj(x) for x in args]
                x, y = sorted((a, b), key=str)
                if T._KIND.get(str(x)) in ('adjoint', 'fundamental'):
                    from .color import KRON
                    return ('s', KRON(x, y))
                return ('s', (T.DELTA if self.p.euclidean else T.G)(x, y))
            if name == 'FSCT':
                from .color import COLT
                i, j, a = [self._slot_obj(x) for x in args]
                return ('s', COLT(a, i, j))
            if name == 'e_':
                objs = [self._slot_obj(x) for x in args]
                e = T._eps_sorted(objs)
                # Minkowski: e_ = -i Eps;  Euclidean: e_ = Eps
                return ('s', e if self.p.euclidean else -I * e)
            if name in ('g_', 'gi_', 'g5_', 'g6_', 'g7_'):
                line = int(args[0])
                if name == 'gi_':
                    return ('g', line, [])
                if name == 'g5_':
                    return ('g', line, [('5',)])
                if name in ('g6_', 'g7_'):
                    return ('g', line, [('16' if name == 'g6_' else '17',)])
                at = []
                for x in args[1:]:
                    if x == '5_':
                        at.append(('5',))
                    elif x in ('6_', '7_'):
                        at.append(('16' if x == '6_' else '17',))
                    else:
                        o = self._obj(x)
                        at.append(('i', o) if T.is_index(o) else ('p', o))
                return ('g', line, at)
            obj = self._obj(name)                       # a vector component p(mu)
            if len(args) == 1 and T.is_momentum(obj):
                return ('s', T.COMP(obj, self._slot_obj(args[0])))
            raise ParseError("unexpected function in FORM output: %s" % s)
        if '.' in s:
            a, b = s.split('.')
            x, y = sorted((self._obj(a), self._obj(b)), key=str)
            return ('s', T.DOT(x, y))
        return ('s', SR(self._obj(s)))

    # -- whole expressions
    def terms(self, body):
        """[(Sage coefficient, {line: [atoms]})] for a printed expression body."""
        if body in ('', '0'):
            return []
        out = []
        fast = not _ARG_OPS.search(body)            # no + - * inside brackets: split with regexes
        for t in (_TERM.split(body) if fast else _split_top(body, '+-')):
            if not t:
                continue
            sign = 1
            if t[0] in '+-':
                sign = -1 if t[0] == '-' else 1
                t = t[1:]
            factors = t.split('*') if fast else _split_top(t, '*')
            sc, lines = [], {}
            coef = Integer(sign)
            for f in factors:
                if f[0].isdigit() and '/' in f:
                    coef *= QQ(f)
                    continue
                a = self.atom(f)
                if a[0] == 's':
                    sc.append(a[1])
                else:
                    lines.setdefault(a[1], []).extend(a[2])
            c = coef * prod(sc) if sc else SR(coef)
            out.append((c, lines))
        return out

    def scalar(self, body):
        ts = self.terms(body)
        if any(l for _, l in ts):
            raise ParseError("Dirac matrices left in a result that should be a number")
        if not ts:
            return SR(0)
        return SR(0).add(*[c for c, _ in ts])


def _depth_ok(s):
    d = 0
    for ch in s:
        d += ch == '('
        d -= ch == ')'
        if d < 0:
            return False
    return d == 0
