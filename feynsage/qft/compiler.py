r"""
Translate feynsage objects into a FORM program (internal).

The walk is over Sage expression trees and DiracExpr terms, never over strings.  Every name gets a
FORM name: the Sage name itself when it is plain ([A-Za-z][A-Za-z0-9]*) and not taken by an object of
another kind, otherwise a generated one.  Pieces of a coefficient that FORM cannot hold as a
polynomial (sqrt(s), log(m), 1/(s - m^2), floats, other functions) are sent as fresh FORM symbols and
put back by the parser.
"""
import re
from sage.all import SR, QQ, Integer
from . import tensors as T

_PLAIN = re.compile(r'[A-Za-z][A-Za-z0-9]*$')


def _colt():
    from .color import COLT
    return COLT(SR.wild(0), SR.wild(1), SR.wild(2))


def _colf():
    from .color import COLF
    return COLF(SR.wild(0), SR.wild(1), SR.wild(2))


class CompileError(TypeError):
    pass


def _opname(e):
    op = e.operator()
    return None if op is None else getattr(op, '__name__', str(op))


class Program:
    """Names and declarations of one FORM run."""

    def __init__(self, dim=4, euclidean=False):
        self.euclidean = bool(euclidean)
        self.has_eps = False
        self.has_color = False
        self.form_of = {}            # (kind, sage name) -> FORM name
        self.sage_of = {}            # FORM name -> (kind, Sage object)
        self.taken = set()
        self.opaque_keys = {}        # str(expr) -> FORM name
        self.dim = dim
        self.dim_name = None
        if isinstance(dim, str):
            dim = 4 if dim.strip() == '4' else SR.var(dim.strip())
            self.dim = dim
        if dim != 4:
            d = SR(dim)
            if not d.is_symbol():
                raise CompileError("dim must be 4 or a symbol such as d or D, got %s" % dim)
            self.dim_name = self.name('symbol', d)

    # -- names
    def name(self, kind, x):
        key = (kind, str(x))
        n = self.form_of.get(key)
        if n is not None:
            return n
        s = str(x)
        if not (_PLAIN.match(s) and s not in self.taken and not s.lower().startswith('fs')):
            k = 1
            while 'fs%s%d' % (kind[0], k) in self.taken:
                k += 1
            s = 'fs%s%d' % (kind[0], k)
        self.taken.add(s)
        self.form_of[key] = s
        self.sage_of[s] = (kind, x)
        return s

    def opaque(self, e):
        k = str(e)
        n = self.opaque_keys.get(k)
        if n is None:
            j = len(self.opaque_keys) + 1
            n = 'fsx%d' % j
            while n in self.taken:
                j += 1
                n = 'fsx%d' % j
            self.taken.add(n)
            self.opaque_keys[k] = n
            self.sage_of[n] = ('opaque', e)
        return n

    def declarations(self):
        groups = {'symbol': [], 'vector': [], 'index': [], 'adjoint': [], 'fundamental': []}
        for (kind, _), n in self.form_of.items():
            groups[kind].append(n)
        groups['symbol'] += [n for n in self.opaque_keys.values()]
        out = []
        if self.has_color or groups['adjoint'] or groups['fundamental']:
            from .color import NC
            self.sage_of['FSNC'] = ('symbol', NC)
            out.append('Symbols FSNC,FSNA;')
            out.append('CFunctions FSCT;')
            out.append('Indices fsi1,fsi2,fsi3,fsi4,fsa1,fsa2;')
            if groups['fundamental']:
                out.append('Indices %s;' % ','.join('%s=FSNC' % n for n in sorted(groups['fundamental'])))
            if groups['adjoint']:
                out.append('Indices %s;' % ','.join('%s=FSNA' % n for n in sorted(groups['adjoint'])))
        if groups['symbol']:
            out.append('Symbols %s;' % ','.join(sorted(groups['symbol'])))
        if self.dim_name:
            out.append('Dimension %s;' % self.dim_name)
        if groups['vector']:
            out.append('Vectors %s;' % ','.join(sorted(groups['vector'])))
        if groups['index']:
            out.append('Indices %s;' % ','.join(sorted(groups['index'])))
        return '\n'.join(out)

    # -- scalars and tensors
    def _slot(self, x):
        """An argument of d_, e_ or a component: a Lorentz index or one momentum."""
        if T.is_index(x):
            return self.name('index', x)
        if T.is_momentum(x):
            return self.name('vector', x)
        raise CompileError("%s is neither a Lorentz index nor a momentum (declare it with lorentz_indices or momenta)" % x)

    def _has_tensor_symbols(self, e):
        return any(T.is_index(v) or T.is_momentum(v) for v in e.variables())

    def number(self, e):
        try:
            q = QQ(e)
            return '(%s)' % q if q.denominator() != 1 or q < 0 else str(q)
        except (TypeError, ValueError):
            pass
        try:
            re_, im_ = QQ(e.real_part()), QQ(e.imag_part())
            return '(%s+(%s)*i_)' % (re_, im_)
        except (TypeError, ValueError):
            return self.opaque(e)

    def scalar(self, e):
        """A Sage expression (scalars and the tensors dot, comp, g, delta, Eps, eps) -> FORM text."""
        e = SR(e)
        op = _opname(e)
        if op is None:
            if e.is_symbol():
                if T.is_momentum(e):
                    raise CompileError("the momentum %s stands alone as a number; use dot(%s, q), comp(%s, mu) or slash(%s)" % (e, e, e, e))
                if T.is_index(e):
                    raise CompileError("the Lorentz index %s stands alone as a number; it must sit in gamma(), metric(), comp() or epsilon()" % e)
                return self.name('symbol', e)
            if e.is_numeric():
                return self.number(e)
            return self.opaque(e)               # pi, euler_gamma, ...
        args = e.operands()
        if op in ('add_vararg', 'add'):
            return '(' + '+'.join(self.scalar(a) for a in args) + ')'
        if op in ('mul_vararg', 'mul'):
            return '*'.join(self.scalar(a) for a in args)
        if op == 'pow':
            b, x = args
            if x.is_integer():
                n = int(x)
                if n >= 0 and (b.has(_colt()) or b.has(_colf())):
                    return '*'.join('(%s)' % self.scalar(b) for _ in range(n)) or '1'   # fresh indices each time
                if n >= 0:
                    return '(%s)^%d' % (self.scalar(b), n)
                if b.is_symbol() and not self._has_tensor_symbols(b):
                    return '%s^(%d)' % (self.name('symbol', b), n)
                if T.TENSOR_NAMES.get(_opname(b)) == 'dot' and all(T.is_momentum(a) for a in b.operands()):
                    return '%s^(%d)' % (self.tensor('dot', b.operands()), n)      # FORM keeps 1/p.k
            return self._opaque_scalar(e)
        if op in T.TENSOR_NAMES:
            return self.tensor(T.TENSOR_NAMES[op], args)
        if op == 'T_color':
            self.has_color = True
            a, i, j = args
            return 'FSCT(%s,%s,%s)' % (self.name('fundamental', i), self.name('fundamental', j), self.name('adjoint', a))
        if op == 'f_color':
            from .color import expand_f
            return '(%s)' % self.scalar(expand_f(e))
        if op == 'Kron':
            a, b = args
            k = 'adjoint' if T._KIND.get(str(a)) == 'adjoint' else 'fundamental'
            self.has_color = True
            return 'd_(%s,%s)' % (self.name(k, a), self.name(k, b))
        if any(T.is_index(v) for v in e.variables()):
            raise CompileError("%s: a Lorentz index is inside the function %s" % (e, op))
        if op == 'gamma' and len(args) == 1 and args[0].is_symbol():
            raise CompileError("gamma(%s) is Euler's Gamma function of the symbol %s here, because %s is not a "
                               "Lorentz index.  Declare it first: %s = lorentz_indices('%s')" % ((args[0],) * 5))
        return self._opaque_scalar(e)

    def _opaque_scalar(self, e):
        """A piece FORM cannot hold (sqrt, 1/(a + b), ...) -> one FORM symbol.  It may contain scalar
        products (it is a number) but no free Lorentz index."""
        if any(T.is_index(v) for v in e.variables()):
            raise CompileError("%s: a Lorentz index cannot sit in a denominator, a root or a function" % e)
        if any(T.is_momentum(v) for v in e.variables()):
            bad = [v for v in e.variables() if T.is_momentum(v)]
            if not self._momenta_only_in_dots(e):
                raise CompileError("%s: the momentum %s must appear inside dot(), comp() or slash()" % (e, bad[0]))
        return self.opaque(e)

    def _momenta_only_in_dots(self, e):
        op = _opname(e)
        if op is None:
            return not T.is_momentum(e)
        if T.TENSOR_NAMES.get(op) in ('dot', 'Eps', 'eps'):
            return True
        return all(self._momenta_only_in_dots(a) for a in e.operands())

    def tensor(self, op, args):
        if op == 'dot':
            if all(T.is_momentum(a) for a in args):
                return '%s.%s' % tuple(self.name('vector', a) for a in args)
            return self.scalar(T.dot(*args))
        if op == 'comp':
            p, mu = args
            if T.is_momentum(p):
                return '%s(%s)' % (self.name('vector', p), self._slot(mu))
            return self.scalar(T.comp(p, mu))
        if op in ('g', 'delta'):
            if (op == 'g') == self.euclidean:
                raise CompileError("the metric (%s, %s) is %s but this calculation is %s; do not mix the conventions"
                                   % (args[0], args[1], 'Minkowski' if op == 'g' else 'Euclidean',
                                      'Euclidean' if self.euclidean else 'Minkowski'))
            return 'd_(%s,%s)' % (self._slot(args[0]), self._slot(args[1]))
        if op == 'Eps':
            if not all(T.is_index(a) or T.is_momentum(a) for a in args):
                return self.scalar(T.epsilon(*args))
            self.has_eps = True
            core = 'e_(%s)' % ','.join(self._slot(a) for a in args)
            # Minkowski: Eps = i e_  (FORM's e_ = -i eps^{mu nu rho sigma});  Euclidean: Eps = e_
            return core if self.euclidean else '(i_*%s)' % core
        if op == 'eps':                       # FORM's own e_, as printed by form.compute
            self.has_eps = True
            return 'e_(%s)' % ','.join(self._slot(a) for a in args)
        raise CompileError(op)

    # -- Dirac matrices
    def factor(self, F, line):
        parts = []
        for a, c in F.items:
            if a[0] == '1':
                m = 'gi_(%d)' % line
            elif a[0] == '5':
                m = 'g5_(%d)' % line
            elif a[0] == 'i':
                m = 'g_(%d,%s)' % (line, self.name('index', a[1]))
            else:
                m = 'g_(%d,%s)' % (line, self.name('vector', a[1]))
            if (c - 1).is_trivial_zero():
                parts.append(m)
            else:
                parts.append('%s*%s' % (self.scalar(c), m))
        if not parts:
            return '0'
        return parts[0] if len(parts) == 1 else '(' + '+'.join(parts) + ')'

    def chain(self, factors, line):
        """A product of factors on one FORM line; runs of plain gamma matrices become one g_(...)."""
        if not factors:
            return 'gi_(%d)' % line
        out, run = [], []

        def flush():
            if run:
                out.append('g_(%d,%s)' % (line, ','.join(run)))
                run.clear()
        for F in factors:
            if len(F.items) == 1 and F.items[0][0][0] in ('i', 'p') and (F.items[0][1] - 1).is_trivial_zero():
                a = F.items[0][0]
                run.append(self.name('index' if a[0] == 'i' else 'vector', a[1]))
            else:
                flush()
                out.append(self.factor(F, line))
        flush()
        return '*'.join(out)

    # -- rules
    def rules(self, rules):
        """{dot(p, q): value, p4: p1 + p2 - p3, m: 0} -> (vector ids, other ids) as FORM statements."""
        vec, other = [], []
        for k, val in (rules or {}).items():
            k = SR(k)
            if T.is_momentum(k):
                parts = T.vector_parts(val, "the right side of a momentum rule")
                rhs = '+'.join('%s*%s' % (self.scalar(c), self.name('vector', p)) for p, c in parts) or '0'
                vec.append('id %s = %s;' % (self.name('vector', k), rhs))
            elif T.TENSOR_NAMES.get(_opname(k)) == 'dot':
                a, b = k.operands()
                if not (T.is_momentum(a) and T.is_momentum(b)):
                    raise CompileError("rule %s: write it for single momenta, e.g. dot(p, k)" % k)
                other.append('id %s.%s = %s;' % (self.name('vector', a), self.name('vector', b), self.scalar(val)))
            elif k.is_symbol() and not T.is_index(k):
                other.append('id %s = %s;' % (self.name('symbol', k), self.scalar(val)))
            else:
                raise CompileError("cannot use %s as a rule; rules are dot(p, q): value, a momentum: a sum of momenta, or a symbol: value" % k)
        return vec, other
