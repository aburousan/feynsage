r"""
Reading the strings a user gives feynsage ("l + p", "s/2", "m^2").

Sage's own SR("...") can look a bare name up in the session: after `l = 1`, SR("l") is the
number 1.  Here every name becomes a symbol of that name, except functions (a name followed by
"(") and the constants pi, I, e, euler_gamma, oo, so the result never depends on the session.
"""
import re
from sage.all import SR, sage_eval

_CONSTANTS = {'pi', 'I', 'e', 'euler_gamma', 'oo', 'infinity'}
_NAME = re.compile(r'[A-Za-z_]\w*')


def sr(text):
    """A string as a symbolic expression with every name a symbol; anything else through SR()."""
    if not isinstance(text, str):
        return SR(text)
    loc = {}
    for m in _NAME.finditer(text):
        n = m.group(0)
        if n in _CONSTANTS or text[m.end():].lstrip().startswith('('):
            continue
        loc[n] = SR.var(n)
    return SR(sage_eval(text, locals=loc))


def swap_i(expr, by):
    """expr with the imaginary unit I replaced by `by` (by = -I conjugates the numbers).

    Sage's expr.subs({I: x}) misses I when it is the numeric coefficient of a single monomial:
    (I*x).subs({I: -I}) gives I*x back, while (I*x + I*y).subs({I: -I}) works.  This walks the
    expression tree and changes every number z to re(z) + by*im(z)."""
    from sage.all import I
    e = SR(expr)
    return _swap_i(e, SR(by)) if e.has(I) else e


def _swap_i(e, by):
    from sage.all import I
    if e.is_numeric():
        z = e.pyobject()
        try:
            re_, im_ = z.real(), z.imag()
        except AttributeError:
            return e
        return e if im_ == 0 else SR(re_) + by * SR(im_)
    op = e.operator()
    if op is None:
        return e
    return op(*[_swap_i(a, by) if a.has(I) else a for a in e.operands()])
