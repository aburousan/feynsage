r"""
LaTeX output for Jupyter: results typeset as equations instead of plain text.

    eq(r"B_{00}(s; m, m)", PVB(1, 0, s, m, m))     # a displayed equation
    loop(...)                                      # LoopResult shows itself in LaTeX
    ibp_reduce(...)                                # so does Reduction

Long expressions are written as a Laurent series in epsilon (1/eps^2, 1/eps, finite part) and
the finite part is broken over several lines.
"""
import re
from sage.all import SR, latex

_GREEK = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda',
          'mu', 'nu', 'xi', 'pi', 'rho', 'sigma', 'tau', 'phi', 'chi', 'psi', 'omega']


def _sym(name):
    return '\\' + name if name in _GREEK else name


def structure(label):
    """'g^{mu nu}' -> g^{\\mu\\nu}, 'p^mu q^nu' -> p^{\\mu} q^{\\nu}, '1' -> 1."""
    if label == '1':
        return '1'
    s = label
    s = re.sub(r'(\w+)\^\{(\w+) (\w+)\}', lambda m: '%s^{%s%s}' % (_sym(m.group(1)), _sym(m.group(2)), _sym(m.group(3))), s)
    s = re.sub(r'(\w+)\^(\w+)', lambda m: '%s^{%s}' % (_sym(m.group(1)), _sym(m.group(2))), s)
    return s


def _clean(t):
    return t.replace(r'\mathit{eps}', r'\epsilon').replace(r'\, ', ' ')


import operator as _op


def _is_add(f):
    o = f.operator()
    return o is not None and 'add' in str(o)


def _is_mul(f):
    o = f.operator()
    return o is not None and 'mul' in str(o)


def _is_func(f):
    """A function application (log, DiscB, LogM, C0, polylog, ...) or a power of one."""
    o = f.operator()
    if o is None or _is_add(f) or _is_mul(f):
        return False
    if o == _op.pow:
        return _is_func(f.operands()[0])
    return True


def grouped(expr):
    """[(coefficient, function part)]: the terms of expr collected by their transcendental factor."""
    e = SR(expr).expand()
    terms = e.operands() if _is_add(e) else [e]
    groups, order = {}, []
    for t in terms:
        fac = t.operands() if _is_mul(t) else [t]
        fpart = SR(1)
        for f in fac:
            if _is_func(f):
                fpart *= f
        key = str(fpart)
        if key not in groups:
            groups[key] = [fpart, SR(0)]
            order.append(key)
        groups[key][1] += t / fpart
    out = []
    for k in order:
        fpart, c = groups[k]
        c = c.factor() if not c.is_numeric() else c
        if not SR(c).is_trivial_zero():
            out.append((c, fpart))
    out.sort(key=lambda cf: str(cf[1]) == '1')        # the rational part last
    return out


def _term(c, f):
    cl = _clean(latex(c))
    if str(f) == '1':
        return cl
    fl = _clean(latex(f))
    if SR(c - 1).is_trivial_zero():
        return fl
    if SR(c + 1).is_trivial_zero():
        return '-' + fl
    if _is_add(SR(c)):
        cl = r'\left(%s\right)' % cl
    return cl + r'\,' + fl


def lines(expr, per_line=3):
    """latex of a sum, broken into rows of per_line terms (aligned, continuation rows start with &)."""
    e = SR(expr)
    ops = e.operands() if e.operator() is not None and 'add' in str(e.operator()) else [e]
    if len(ops) <= per_line:
        return _clean(latex(e))
    rows = []
    for i in range(0, len(ops), per_line):
        chunk = ' + '.join(_clean(latex(t)) for t in ops[i:i + per_line])
        rows.append(chunk)
    body = r' \\ & \quad + '.join(rows)
    return body.replace('+ -', '- ')


def laurent(expr, per_line=2):
    """c_{-2}/eps^2 + c_{-1}/eps + c_0: poles first, then the finite part collected by functions
    (one log, DiscB, C0, ... per term, coefficients factored), per_line terms on a row."""
    from .pv import eps
    e = SR(expr).expand()
    pieces = []
    for k, head in ((-2, r'\frac{1}{\epsilon^2}'), (-1, r'\frac{1}{\epsilon}')):
        c = e.coefficient(eps, k) if e.has(eps) else SR(0)
        if not SR(c).is_trivial_zero():
            if SR(c).is_numeric():                     # a number: write c/eps^k as one fraction
                pieces.append(_clean(latex(SR(c) * eps ** k)))
            else:
                pieces.append(r'%s\left(%s\right)' % (head, _clean(latex(SR(c).factor()))))
    c0 = e.coefficient(eps, 0) if e.has(eps) else e
    pieces += [_term(c, f) for c, f in grouped(c0)]
    if not pieces:
        return '0'
    rows = [' + '.join(pieces[i:i + per_line]) for i in range(0, len(pieces), per_line)]
    return r' \\ & \quad + '.join(rows).replace('+ -', '- ')


def eq(lhs, expr, per_line=2):
    """Display  lhs = expr  as an equation in Jupyter (plain text elsewhere)."""
    body = r'\begin{aligned} %s &= %s \end{aligned}' % (lhs, laurent(expr, per_line))
    try:
        from IPython.display import Math, display
        display(Math(body))
    except ImportError:
        print('%s = %s' % (lhs, expr))
