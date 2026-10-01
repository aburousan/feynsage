r"""
A small bridge to FORM (J. Vermaseren) for Dirac algebra.

    trace(['mu','nu','rho','sigma'])          -> Tr(g_mu g_nu g_rho g_sigma) in D dims
    dirac_trace(factors, vectors)             -> one fermion line with slashed momenta and gamma_5
    compute(expr, vectors, lines, rules)      -> traces and contractions written in Sage notation,
                                                 converted to FORM and run (no FORM code to write)
    run_form(code)                            -> raw FORM output, for FORM programs written by hand

The result is returned as a Sage symbolic expression in d_(mu,nu) -> g(mu,nu) and
p.q -> dot(p,q).
"""
import os, re, subprocess, tempfile, shutil
from sage.all import SR, function

def _find_form():
    """FEYNSAGE_FORM, then the PATH, then the places install.sh and package managers use."""
    cands = [os.environ.get('FEYNSAGE_FORM'), shutil.which('form'),
             os.path.expanduser('~/.local/bin/form'), os.path.expanduser('~/bin/form'),
             '/opt/homebrew/bin/form', '/usr/local/bin/form', '/usr/bin/form']
    for c in cands:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return 'form'


FORM = _find_form()


def run_form(code):
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, 'job.frm')
        open(f, 'w').write(code)
        r = subprocess.run([FORM, '-q', f], capture_output=True, text=True, cwd=tmp)
        if r.returncode != 0:
            raise RuntimeError(r.stdout + r.stderr)
        return r.stdout


def trace(indices, dim='D'):
    """Trace of a product of gamma matrices gamma_{i1} ... gamma_{in} in dim dimensions."""
    idx = ','.join(indices)
    code = f"""
Symbol {dim};
Dimension {dim};
Indices {idx};
Local T = g_(1,{idx});
tracen,1;
Print +s;
.end
"""
    out = run_form(code)
    body = out.split('T =', 1)[1].split(';', 1)[0]
    body = re.sub(r'\s+', '', body)
    body = re.sub(r'd_\(([^,]+),([^)]+)\)', r'g(\1,\2)', body)
    g = function('g')
    names = {i: SR.var(i) for i in indices}
    names['g'] = g
    names[dim] = SR.var(dim)
    return SR(eval(body.replace('^', '**'), {}, names))


def dirac_trace(factors, vectors, dim=4, levi_civita="form"):
    r"""
    Trace of a product of Dirac matrices, done by FORM.

    Each factor is a string:
      'g5'                  gamma_5
      an index name 'mu'    gamma^mu (any name that is not a vector)
      a linear expression   slash of the vector part plus the rest times 1,
                            e.g. 'q - x*k1 + y*k2 + m'  ->  q/ - x k1/ + y k2/ + m
    vectors: the names that are momenta.  dim = 4 uses FORM's trace4 (needed with g5),
    any other value is a symbolic dimension (no g5 then).

    levi_civita="usual" writes the answer with the usual epsilon of Peskin and Schroeder,
    Eps (= i e_, eps^0123 = -1); the pion trace is then 4*I*m*Eps(k1, k2, mu, nu).

    Returns a Sage expression with g(mu, nu) for the metric, dot(p, q) for p.q,
    comp(p, mu) for the component p^mu and eps(a, b, c, d) for FORM's e_.  FORM's e_ is -i times the usual Levi-Civita tensor
    epsilon (the one of FeynCalc's Eps): checked on the pion trace, where FORM gives
    4 m e_(k1,k2,mu,nu) and FeynCalc -4 i m Eps(mu,nu,k1,k2).

        dirac_trace(['g5', 'q - x*k1 + y*k2 + m', 'mu', 'q + (1-x)*k1 + y*k2 + m', 'nu',
                     'q - x*k1 - (1-y)*k2 + m'], vectors=['q', 'k1', 'k2'])
        -> 4*m*eps(k1, k2, mu, nu)
    """
    from sage.all import SR, var
    vecs = [SR.var(v) for v in vectors]
    indices, symbols, fac = [], set(), []
    for f in factors:
        f = f.strip()
        if f == 'g5':
            fac.append('g_(1,5_)')
            continue
        if re.fullmatch(r'[A-Za-z]\w*', f) and f not in vectors:
            indices.append(f)
            fac.append('g_(1,%s)' % f)
            continue
        e = SR(f).expand()
        parts = []
        rest = e
        for v in vecs:
            c = e.coefficient(v)
            rest = rest - c * v
            if c != 0:
                parts.append('(%s)*g_(1,%s)' % (c, v))
                symbols.update(str(s) for s in c.variables())
        rest = rest.expand()
        if rest != 0:
            parts.append('(%s)*gi_(1)' % rest)
            symbols.update(str(s) for s in rest.variables())
        fac.append('(' + ' + '.join(parts) + ')' if parts else '0')
    decl = []
    if dim != 4:
        decl += ['Symbol %s;' % dim, 'Dimension %s;' % dim]
        symbols.discard(str(dim))
    if symbols:
        decl.append('Symbols %s;' % ','.join(sorted(symbols)))
    if vectors:
        decl.append('Vectors %s;' % ','.join(vectors))
    if indices:
        decl.append('Indices %s;' % ','.join(sorted(set(indices))))
    code = '\n'.join(decl) + '\nLocal T = ' + '*'.join(fac) + ';\n' + \
           ('trace4,1;\n' if dim == 4 else 'tracen,1;\n') + 'Format nospaces;\nPrint +s;\n.end\n'
    out = run_form(code)
    body = re.split(r'\bT\s*=', out, maxsplit=1)[1].split(';', 1)[0]
    body = re.sub(r'\s+', '', body)
    body = re.sub(r'd_\(([^,]+),([^)]+)\)', r'g(\1,\2)', body)
    body = body.replace('e_(', 'eps(')
    body = re.sub(r'\b(\w+)\.(\w+)\b', r'dot(\1,\2)', body)
    for v in vectors:                                   # p(mu) is the component p^mu
        body = re.sub(r'(?<![\w,(])%s\(' % v, 'comp(%s,' % v, body)
        body = re.sub(r'(?<=[(,*+-])%s\((?=\w+\))' % v, 'comp(%s,' % v, body)
    names = {'g': function('g'), 'eps': function('eps'), 'dot': function('dot'), 'comp': function('comp')}
    for n in set(re.findall(r'[A-Za-z_]\w*', body)) - set(names):
        names[n] = SR.var(n)
    return _levi_civita(SR(eval(body.replace('^', '**'), {}, names)), levi_civita)


EPS = function('Eps', latex_name=r'\varepsilon')        # the usual Levi-Civita symbol


def _levi_civita(expr, convention):
    """FORM's e_ (printed eps) -> the usual epsilon Eps = i e_, if convention == "usual"."""
    if convention == "form":
        return expr
    if convention != "usual":
        raise ValueError('levi_civita must be "form" or "usual"')
    I_ = SR(1).parent()('I')
    return expr.substitute_function(function('eps'), lambda *a: -I_ * EPS(*a)).expand()


def _line(factors, vectors, n, symbols, indices):
    """One fermion line as FORM: 'g5' -> g_(n,5_), an index 'mu' -> g_(n,mu), a linear
    expression 'q - x*k1 + m' -> (g_(n,q) - x*g_(n,k1) + m*gi_(n))."""
    vecs = [SR.var(v) for v in vectors]
    out = []
    for f in factors:
        f = f.strip()
        if f == 'g5':
            out.append('g_(%d,5_)' % n)
            continue
        if re.fullmatch(r'[A-Za-z]\w*', f) and f not in vectors:
            indices.add(f)
            out.append('g_(%d,%s)' % (n, f))
            continue
        e = SR(f).expand()
        parts, rest = [], e
        for v in vecs:
            c = e.coefficient(v)
            rest = rest - c * v
            if c != 0:
                parts.append('(%s)*g_(%d,%s)' % (c, n, v))
                symbols.update(str(x) for x in c.variables())
        rest = rest.expand()
        if rest != 0:
            parts.append('(%s)*gi_(%d)' % (rest, n))
            symbols.update(str(x) for x in rest.variables())
        out.append('(' + ' + '.join(parts) + ')' if parts else '0')
    return '*'.join(out)


def _sage_to_form(expr, vectors, symbols, indices):
    """A Sage-style expression -> FORM: eps(...) -> e_(...), g(a,b) -> d_(a,b), dot(p,q) and
    p.q -> p.q, comp(p,mu) and p(mu) -> p(mu); collects the indices and symbols it meets."""
    t = str(expr).replace('**', '^')
    t = re.sub(r'\bdot\(\s*(\w+)\s*,\s*(\w+)\s*\)', r'\1.\2', t)
    t = re.sub(r'\bcomp\(\s*(\w+)\s*,\s*(\w+)\s*\)', r'\1(\2)', t)
    t = re.sub(r'\beps\(', 'e_(', t)
    t = re.sub(r'\bg\(', 'd_(', t)
    for args in re.findall(r'\b(?:e_|d_)\(([^)]*)\)', t):
        for a in args.split(','):
            a = a.strip()
            if a and a not in vectors and not re.fullmatch(r'-?\d+', a):
                indices.add(a)
    for v in vectors:
        for a in re.findall(r'\b%s\(\s*(\w+)\s*\)' % re.escape(v), t):
            if a not in vectors:
                indices.add(a)
    for name in re.findall(r'\b[A-Za-z]\w*\b', re.sub(r'\b\w+\.\w+\b', ' ', t)):
        if name not in vectors and name not in indices and name not in ('e_', 'd_', 'i_'):
            symbols.add(name)
    return t


def _parse(out, name, vectors):
    """FORM's printed expression `name` -> Sage, with g, eps, dot and comp as in dirac_trace."""
    body = re.split(r'\b%s\s*=' % name, out, maxsplit=1)[1].split(';', 1)[0]
    body = re.sub(r'\s+', '', body)
    body = re.sub(r'd_\(([^,]+),([^)]+)\)', r'g(\1,\2)', body)
    body = body.replace('e_(', 'eps(').replace('i_', 'I')
    body = re.sub(r'\b(\w+)\.(\w+)\b', r'dot(\1,\2)', body)
    for v in vectors:
        body = re.sub(r'(?<![\w,(])%s\(' % v, 'comp(%s,' % v, body)
        body = re.sub(r'(?<=[(,*+-])%s\((?=\w+\))' % v, 'comp(%s,' % v, body)
    names = {'g': function('g'), 'eps': function('eps'), 'dot': function('dot'), 'comp': function('comp'), 'I': SR(1).parent()('I')}
    for n_ in set(re.findall(r'[A-Za-z_]\w*', body)) - set(names):
        names[n_] = SR.var(n_)
    return SR(eval(body.replace('^', '**'), {}, names))


def compute(expr="1", vectors=(), lines=(), rules=None, dim=4, show_code=False, levi_civita="form"):
    r"""
    Let FORM do Dirac traces and index contractions without writing FORM code.

    expr     a product in Sage-like notation: eps(mu,nu,rho,sigma) (FORM's Levi-Civita e_, which
             is -i times the usual epsilon), g(mu,nu) (the metric), dot(p,q) or p.q, components
             p(mu) or comp(p,mu), and ordinary symbols
    vectors  the names that are momenta
    lines    fermion lines to trace, each a list of factors as in dirac_trace:
             'g5', an index 'mu', or a linear expression 'p + m' (slashed momentum plus mass)
    rules    scalar products to put in after contracting, e.g. {"k1.k1": 0, "k1.k2": "s/2"}
    dim      4 (trace4, needed with g5) or a symbol such as 'D' (tracen)
    show_code=True prints the FORM program that was run.
    levi_civita="usual": eps(...) in expr and Eps(...) in the answer are the usual Levi-Civita
             symbol of Peskin and Schroeder (eps^0123 = -1), so that Tr[g5 g^mu g^nu g^rho g^sigma]
             = -4 i Eps(mu,nu,rho,sigma) and Eps.Eps = -24; the default "form" keeps FORM's
             e_ = -i epsilon, printed eps.

    Every repeated index is summed (contract).  Returns a Sage expression with g, eps, dot, comp.

        compute("eps(mu,nu,rho,sigma)*eps(mu,nu,al,be)*k1(rho)*k2(sigma)*k1(al)*k2(be)",
                vectors=["k1", "k2"], rules={"k1.k1": 0, "k2.k2": 0})        # -> -2*dot(k1,k2)^2
        compute(lines=[["pp", "mu", "p", "nu"], ["k", "mu", "kp", "nu"]],
                vectors=["p", "pp", "k", "kp"])                               # e+ e- -> mu+ mu-
    """
    vectors = list(vectors)
    symbols, indices = set(), set()
    parts = [] if str(expr).strip() in ("1", "") else [_sage_to_form(expr, vectors, symbols, indices)]
    for n, fac in enumerate(lines, start=1):
        parts.append(_line(fac, vectors, n, symbols, indices))
    ids = []
    for k, v in (rules or {}).items():
        k = re.sub(r'\bdot\(\s*(\w+)\s*,\s*(\w+)\s*\)', r'\1.\2', str(k))
        v = _sage_to_form(v, vectors, symbols, indices)
        ids.append('id %s = %s;' % (k, v))
    if dim != 4:
        symbols.discard(str(dim))
    symbols -= set(vectors) | indices
    decl = []
    if dim != 4:
        decl += ['Symbol %s;' % dim, 'Dimension %s;' % dim]
    if symbols:
        decl.append('Symbols %s;' % ','.join(sorted(symbols)))
    if vectors:
        decl.append('Vectors %s;' % ','.join(vectors))
    if indices:
        decl.append('Indices %s;' % ','.join(sorted(indices)))
    code = '\n'.join(decl) + '\nLocal E = ' + ('*'.join(parts) if parts else '1') + ';\n'
    for n in range(1, len(lines) + 1):
        code += ('trace4,%d;\n' % n) if dim == 4 else ('tracen,%d;\n' % n)
    code += 'contract;\n' + ''.join(i + '\n' for i in ids) + 'Format nospaces;\nPrint +s;\n.end\n'
    if show_code:
        print(code)
    res = _parse(run_form(code), 'E', vectors)
    if levi_civita == "usual":                 # usual epsilon = i e_ in the input
        res = res * SR(1).parent()('I')**len(re.findall(r'\beps\(', str(expr)))
    return _levi_civita(res, levi_civita)


def to_loop(expr):
    """Turn a dirac_trace result into a numerator string for pv.loop():
    comp(l, mu) -> l^mu, dot(a, b) -> a.b, g(mu, nu) -> g(mu,nu)."""
    t = str(expr)
    t = re.sub(r'comp\((\w+),\s*(\w+)\)', r'\1^\2', t)
    t = re.sub(r'dot\((\w+),\s*(\w+)\)', r'\1.\2', t)
    return t
