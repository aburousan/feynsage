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
