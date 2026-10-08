r"""
Read the reduction tables that Kira writes (kira_<family>.m, Mathematica syntax):

    table = read_kira("kira_targets.m")      # {(2, 1): expression in the masters, ...}

Every integral fam[a1,...] becomes the Sage symbol fam(a1,...), the same as Family.integral, so the
table can be compared with ibp_reduce directly.
"""
import re
from sage.all import SR
from ._parse import sr

from sage.repl.preparse import preparse as _preparse
from sage.all import Integer as _Integer, RealNumber as _RealNumber
_NUMS = {'Integer': _Integer, 'RealNumber': _RealNumber}    # exact rationals: 8/9 is not a float



def _top_level_split(text):
    """The entries of a Mathematica list: split at the commas that are outside every bracket."""
    out, depth, start = [], 0, 0
    for i, ch in enumerate(text):
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        elif ch == ',' and depth == 0:
            out.append(text[start:i]); start = i + 1
    out.append(text[start:])
    return out


def read_kira(path, family_name=None):
    text = open(path).read()
    text = text.strip().lstrip('{').rstrip('}')
    out = {}
    for entry in _top_level_split(text):
        if '->' not in entry:
            continue
        lhs, rhs = entry.split('->', 1)
        lhs = lhs.strip()
        m = re.match(r'(\w+)\[([^\]]*)\]', lhs)
        name = family_name or m.group(1)
        idx = tuple(int(v) for v in m.group(2).split(','))
        rhs = re.sub(r'(\w+)\[([^\]]*)\]', lambda mm: '%s(%s)' % (family_name or mm.group(1), mm.group(2)), rhs)
        rhs = ' '.join(rhs.split())
        from sage.symbolic.function_factory import function as _fn
        loc = {name: _fn(name)}
        for n in set(re.findall(r'[A-Za-z_]\w*', rhs)) - {name}:
            loc[n] = SR.var(n)
        out[idx] = SR(eval(_preparse(rhs.replace('^', '**')), _NUMS, loc))
    return out
