r"""
A small bridge to FORM (J. Vermaseren) for Dirac algebra.

    trace(['mu','nu','rho','sigma'])          -> Tr(g_mu g_nu g_rho g_sigma) in D dims
    run_form(code)                            -> raw FORM output

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


def dirac_trace(factors, vectors, dim=4):
    r"""
    Trace of a product of Dirac matrices, done by FORM.

    Each factor is a string:
      'g5'                  gamma_5
      an index name 'mu'    gamma^mu (any name that is not a vector)
      a linear expression   slash of the vector part plus the rest times 1,
                            e.g. 'q - x*k1 + y*k2 + m'  ->  q/ - x k1/ + y k2/ + m
    vectors: the names that are momenta.  dim = 4 uses FORM's trace4 (needed with g5),
    any other value is a symbolic dimension (no g5 then).

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
    return SR(eval(body.replace('^', '**'), {}, names))


def to_loop(expr):
    """Turn a dirac_trace result into a numerator string for pv.loop():
    comp(l, mu) -> l^mu, dot(a, b) -> a.b, g(mu, nu) -> g(mu,nu)."""
    t = str(expr)
    t = re.sub(r'comp\((\w+),\s*(\w+)\)', r'\1^\2', t)
    t = re.sub(r'dot\((\w+),\s*(\w+)\)', r'\1.\2', t)
    return t
