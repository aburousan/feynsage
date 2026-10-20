# Reading FeynCalc's InputForm text in the tests.  A large sum is cut at its top-level + and - signs and evaluated
# a few hundred terms at a time: Python compiles a+b+c+... as one deeply nested expression and older versions
# (3.11) stop with a RecursionError on the largest formulas.
from sage.all import sage_eval


def terms(txt):
    """The top-level terms of a sum, signs kept."""
    out, depth, start = [], 0, 0
    for k, ch in enumerate(txt):
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        elif ch in '+-' and depth == 0 and k > start and txt[k - 1] not in '*/^(eE':
            out.append(txt[start:k])
            start = k
    out.append(txt[start:])
    return [t for t in out if t.strip()]


def mathematica(txt):
    """InputForm to Python syntax (brackets of functions and lists, powers)."""
    return txt.replace('Sqrt[', 'sqrt(').replace('[', '(').replace(']', ')').replace('{', '[').replace('}', ']').replace('^', '**')


def evaluate(txt, names, chunk=300):
    """sum of the terms of txt (Python syntax) evaluated with names, chunk terms at a time."""
    ts = terms(txt.strip())
    total = 0
    for i in range(0, len(ts), chunk):
        total += sage_eval(' '.join(ts[i:i + chunk]), locals=names)
    return total


def items(txt):
    """The elements of a top-level list {a, b, ...}."""
    txt = txt.strip()[1:-1]
    out, depth, start = [], 0, 0
    for k, ch in enumerate(txt):
        depth += ch in '([{'
        depth -= ch in ')]}'
        if ch == ',' and depth == 0:
            out.append(txt[start:k])
            start = k + 1
    return out + [txt[start:]]
