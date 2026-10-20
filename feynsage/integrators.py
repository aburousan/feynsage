r"""
Integrals with a check.  Sage's integrate() hands the work to one backend (Maxima unless told otherwise) and
returns whatever comes back: a wrong branch of a logarithm, an answer that silently assumed the sign of a
parameter, or the integral itself unevaluated.  integrate_checked() tries several backends in turn and
returns the first answer that agrees with a numerical integral at random values of the other symbols.  If
none agrees it raises IntegrationError: a wrong formula is never returned.

    integrate_checked(f, x, a, b)       Int_a^b f dx (a, b may be +-oo); without a, b the antiderivative,
                                        checked by differentiating it at random points
    set_integrators(names)              the order in which the backends are tried
    available_integrators()             the backends this installation has

Backends, in the default order (the ones that are not installed are skipped):

    "polynomial"  exact, no backend: a polynomial in x is integrated term by term
    "wolfram"     Mathematica or the free Wolfram Engine, if a kernel is found (FEYNSAGE_WOLFRAM = path to
                  WolframKernel, or `math`, `wolfram`, `WolframKernel` on the PATH, or `wolframscript`)
    "fricas"      FriCAS (`sage -i fricas`), the strongest free integrator (Risch algorithm)
    "maxima"      ships with Sage (Sage's default)
    "sympy"       ships with Sage
    "giac"        ships with Sage; an answer that still contains the integral is not accepted
    "antiderivative"  for a definite integral the free backends cannot do (typically one that gives
                  dilogarithms): the antiderivative from Maxima or SymPy, then one-sided limits at the end
                  points, taken after expand_log so that singular logarithms cancel

Each backend has FEYNSAGE_INTEGRATE_TIMEOUT seconds (default 120); one that takes longer is stopped and the
next is tried.  The environment variable FEYNSAGE_INTEGRATORS ("maxima,sympy", ...) sets the order for a
whole session.
"""
import os
import random
import shutil
import subprocess
import warnings

from sage.all import SR, QQ, CC, Infinity, integrate, assume, forget, assumptions

_ALL = ("polynomial", "wolfram", "fricas", "maxima", "sympy", "giac", "antiderivative")
_order = None
_wolfram = None          # a running _Kernel, or the string 'wolframscript', or False


class IntegrationError(ArithmeticError):
    """No backend gave an answer that passes the numerical check."""


def set_integrators(names):
    """Try the backends in this order, e.g. set_integrators(["wolfram", "maxima"]).  None restores the default."""
    global _order
    if names is None:
        _order = None
        return
    names = [n.strip().lower() for n in ([names] if isinstance(names, str) else names)]
    bad = [n for n in names if n not in _ALL]
    if bad:
        raise ValueError("unknown integrator(s) %s; known: %s" % (bad, ", ".join(_ALL)))
    _order = names


def _wanted():
    if _order is not None:
        return list(_order)
    env = os.environ.get("FEYNSAGE_INTEGRATORS")
    if env:
        set_integrators(env.split(","))
        return list(_order)
    return list(_ALL)


def _kernel_command():
    """The command that starts a Wolfram kernel, or None."""
    path = os.environ.get("FEYNSAGE_WOLFRAM")
    if path:
        return path
    for name in ("math", "wolfram", "WolframKernel"):
        p = shutil.which(name)
        if p:
            return p
    for app in ("/Applications/Wolfram.app/Contents/MacOS/WolframKernel",      # macOS, not on the PATH
                "/Applications/Mathematica.app/Contents/MacOS/WolframKernel",
                "/Applications/Wolfram Engine.app/Contents/MacOS/WolframKernel"):
        if os.path.exists(app):
            return app
    return None


class _Kernel:
    """A Wolfram kernel kept running on a pipe; each command is answered up to an end marker."""
    END = "FEYNSAGE-END"

    def __init__(self, cmd):
        self.p = subprocess.Popen([cmd, "-noprompt"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, bufsize=1)

    def eval(self, code):
        self.p.stdin.write('WriteString[$Output, ToString[%s] <> "\\n%s\\n"]\n' % (code.replace("\n", " "), self.END))
        self.p.stdin.flush()
        lines = []
        while True:
            line = self.p.stdout.readline()
            if not line:
                raise RuntimeError("the Wolfram kernel stopped")
            if line.strip() == self.END:
                return "".join(lines).strip()
            lines.append(line)


def _wolfram_session():
    """Start (once) the Wolfram kernel; 'wolframscript' if only that is there; False if there is none."""
    global _wolfram
    if _wolfram is None:
        cmd = _kernel_command()
        if cmd:
            try:
                K = _Kernel(cmd)
                if K.eval("1 + 1") == "2":
                    _wolfram = K
            except Exception:
                _wolfram = None
        if _wolfram is None:
            _wolfram = "wolframscript" if shutil.which("wolframscript") else False
    return _wolfram


def _has(name):
    if name in ("polynomial", "maxima", "sympy", "giac", "antiderivative"):
        return True
    if name == "fricas":
        try:
            from sage.features.fricas import FriCAS
            return bool(FriCAS().is_present())
        except Exception:
            return False
    if name == "wolfram":
        return bool(_kernel_command() or shutil.which("wolframscript"))
    return False


def available_integrators():
    """The backends that are installed, in the order they are tried."""
    return [n for n in _wanted() if _has(n)]


# ---------------------------------------------------------------------------- the backends
def _polynomial(f, x, a, b):
    if not f.is_polynomial(x):
        return None
    cf = f.expand().coefficients(x)
    if a is None:
        return sum(c * x**(n + 1) / (n + 1) for c, n in cf)
    if any(SR(e).is_infinity() or SR(e).has(Infinity) for e in (a, b)):
        return None
    return sum(c * (SR(b)**(n + 1) - SR(a)**(n + 1)) / (n + 1) for c, n in cf)


def _sage_backend(name):
    def run(f, x, a, b):
        r = integrate(f, x, algorithm=name) if a is None else integrate(f, x, a, b, algorithm=name)
        return None if "integrate(" in str(r) or "limit(" in str(r) else r
    return run


def _antiderivative(f, x, a, b):
    """Int_a^b f dx = F(b-) - F(a+) with F from Maxima or SymPy.  At a finite end point e the limit is taken
    as x = e -+ h, h -> 0+, after expand_log: Maxima then sees the logarithms that cancel."""
    from sage.all import limit
    if a is None:
        return None
    h = SR.var("fs_end_h")
    for name in ("maxima", "sympy"):
        try:
            F = integrate(f, x, algorithm=name)
        except Exception:
            continue
        if "integrate(" in str(F):
            continue
        ends = []
        for e, sign in ((SR(b), -1), (SR(a), 1)):
            tries = ([lambda: limit(F.expand_log(), **{str(x): e})] if e.is_infinity() else
                     [lambda: limit(F.subs({x: e + sign*h}).expand_log(), fs_end_h=0, dir="plus"),
                      lambda: limit(F.expand_log(), dir="minus" if sign < 0 else "plus", **{str(x): e})])
            v = None
            for t in tries:
                try:
                    v = t()
                except Exception:
                    v = None
                if v is not None and "limit(" not in str(v) and not any(
                        w in str(v) for w in ("und", "ind", "infinity")):
                    break
                v = None
            if v is None:
                break
            ends.append(v)
        if len(ends) == 2:
            return ends[0] - ends[1]
    return None


def _wolfram_backend(f, x, a, b):
    from sage.interfaces.mathematica import symbolic_expression_from_mathematica_string as from_m
    S = _wolfram_session()
    if not S:
        return None
    rng = x._mathematica_init_() if a is None else "{%s, %s, %s}" % (
        x._mathematica_init_(), SR(a)._mathematica_init_(), SR(b)._mathematica_init_())
    mine = set(f.variables()) | set(SR(a).variables() if a is not None else ()) | set(SR(b).variables() if b is not None else ())
    facts = []
    for c in assumptions():                    # only the relations among the symbols of this integral
        try:
            c = SR(c)
            if c.is_relational() and set(c.variables()) <= mine:
                facts.append(c._mathematica_init_())
        except (TypeError, ValueError):
            pass
    opt = ", Assumptions -> %s" % " && ".join(facts) if facts else ""
    code = "ToString[Integrate[%s, %s%s], InputForm]" % (f._mathematica_init_(), rng, opt)
    if S == "wolframscript":
        out = subprocess.run(["wolframscript", "-code", code], capture_output=True, text=True, timeout=600).stdout
    else:
        out = S.eval(code)
    out = str(out).strip()
    if not out or "Integrate[" in out or "ConditionalExpression" in out or "$Failed" in out:
        return None
    return from_m(out)


_RUN = {"polynomial": _polynomial, "wolfram": _wolfram_backend, "fricas": _sage_backend("fricas"),
        "maxima": _sage_backend("maxima"), "sympy": _sage_backend("sympy"), "giac": _sage_backend("giac"),
        "antiderivative": _antiderivative}


def _run(name, f, x, a, b):
    """One backend, stopped after FEYNSAGE_INTEGRATE_TIMEOUT seconds (a stopped Wolfram kernel is restarted)."""
    global _wolfram
    from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
    try:
        alarm(float(os.environ.get("FEYNSAGE_INTEGRATE_TIMEOUT", 120)))
        return _RUN[name](f, x, a, b)
    except AlarmInterrupt:
        if name == "wolfram" and isinstance(_wolfram, _Kernel):
            _wolfram.p.kill()
            _wolfram = None
        raise TimeoutError("%s took too long" % name)
    finally:
        cancel_alarm()


# ---------------------------------------------------------------------------- the check
def _holds(facts, point):
    for c in facts:
        try:
            if not bool(SR(c).subs(point)):
                return False
        except Exception:
            return False
    return True


def _points(syms, facts, n, rnd):
    """n random rational points for the symbols that obey the assumptions on them.  Point k is drawn at
    the scale k mod 3 (about 0.5, about 2, 10 to 100): an answer that is right only for parameters of order
    one (a lost branch, or floating-point coefficients that cancel badly at large values) is caught."""
    ranges = [(QQ(1) / 7, QQ(6) / 7), (QQ(1), QQ(3)), (QQ(10), QQ(100)), (QQ(-3), QQ(-1) / 3), (QQ(-100), QQ(-10))]
    mine = []
    for c in facts:
        try:
            if set(SR(c).variables()) <= set(syms):
                mine.append(c)
        except (TypeError, ValueError):                  # "n is integer" and the like
            pass
    pts = []
    for k in range(n):
        order = [k % 3] + [i for i in range(len(ranges)) if i != k % 3]
        for tries in range(200):
            lo, hi = ranges[order[min(tries // 20, len(order) - 1)]]
            p = {s: lo + (hi - lo) * QQ(rnd.randint(1, 997)) / 998 for s in syms}
            if _holds(mine, p):
                pts.append(p)
                break
    return pts


def _numeric_definite(f, x, a, b):
    """The integral by GSL's adaptive QAGS (end points are never evaluated, so log(x) at 0 is fine), real and
    imaginary parts apart; mpmath at 30 digits when GSL is not sure of 12 digits."""
    import math
    import mpmath
    from sage.all import numerical_integral, CDF, ComplexField
    from sage.ext.fast_callable import fast_callable
    lo, hi = SR(a), SR(b)
    lo = -Infinity if lo == -Infinity else float(lo.n())
    hi = Infinity if hi == Infinity else float(hi.n())
    g = fast_callable(SR(f), vars=[x], domain=CDF)
    out, sure = [], True
    for part in (lambda t: g(t).real(), lambda t: g(t).imag()):
        v, err = numerical_integral(part, lo, hi, eps_abs=1e-14, eps_rel=1e-12, max_points=200)
        out.append(v)
        sure = sure and math.isfinite(v) and err <= 1e-11 * max(1.0, abs(v))
    if sure:
        return complex(out[0], out[1])
    F = ComplexField(120)
    h = fast_callable(SR(f), vars=[x], domain=F)
    with mpmath.workdps(30):
        lo_ = -mpmath.inf if lo == -Infinity else mpmath.mpf(lo)
        hi_ = mpmath.inf if hi == Infinity else mpmath.mpf(hi)
        v = mpmath.quad(lambda t: mpmath.mpc(*[mpmath.mpf(str(c)) for c in
                                               (lambda z: (z.real(), z.imag()))(h(F(str(t))))]), [lo_, hi_])
    return complex(v)


def _close(u, v, tol):
    return abs(u - v) <= tol * max(1.0, abs(u), abs(v))


def _passes(F, f, x, a, b, points, tol):
    for p in points:
        try:
            if a is None:                                    # antiderivative: F' = f at two points
                dF = SR(F).diff(x).subs(p)
                for t in (QQ(2) / 7, QQ(5) / 3):
                    u, v = complex(CC(dF.subs({x: t}).n(80))), complex(CC(SR(f).subs(p).subs({x: t}).n(80)))
                    if not _close(u, v, tol):
                        return False
            else:
                u = complex(CC(SR(F).subs(p).n(80)))
                v = _numeric_definite(SR(f).subs(p), x, SR(a).subs(p), SR(b).subs(p))
                if not _close(u, v, tol):
                    return False
        except (ValueError, TypeError, ZeroDivisionError, ArithmeticError):
            return False
    return True


def integrate_checked(f, x, a=None, b=None, assume_=(), check=True, points=3, tol=1e-8, seed=1,
                      info=False, integrators=None):
    """
    Int_a^b f dx (or the antiderivative when a, b are None), tried with each backend in turn (see the module
    doc) and checked numerically: at `points` random rational values of the other symbols (chosen to obey
    the assumptions in force and those in assume_) the answer must agree with mpmath's quadrature to `tol`
    (an antiderivative must differentiate back to f).  The first answer that passes is returned; with
    info=True as (answer, backend).  Raises IntegrationError if no backend passes.  check=False returns the
    first answer that is not an unevaluated integral, unchecked.
    """
    f, x = SR(f), SR(x)
    extra = [SR(c) for c in assume_]
    if extra:
        assume(*extra)
    try:
        facts = list(assumptions())
        syms = sorted((s for s in f.variables() if s != x), key=str)
        for e in (a, b):
            if e is not None:
                syms = sorted(set(syms) | set(SR(e).variables()), key=str)
        pts = _points(syms, facts, points, random.Random(seed)) if syms else [{}]
        if check and syms and not pts:
            warnings.warn("integrate_checked: no point obeys the assumptions; the answer is not checked")
        names = integrators or _wanted()
        tried = []
        for name in names:
            if not _has(name):
                continue
            try:
                F = _run(name, f, x, a, b)
            except Exception as err:                       # a backend that fails or asks a question
                tried.append("%s (%s)" % (name, type(err).__name__))
                continue
            if F is None:
                tried.append("%s (no answer)" % name)
                continue
            if not check or not pts or _passes(F, f, x, a, b, pts, tol):
                return (F, name) if info else F
            tried.append("%s (failed the numerical check)" % name)
        raise IntegrationError("no integrator gave a checked answer for Int %s d%s%s; tried: %s.  Install a "
                               "stronger one (FriCAS: `sage -i fricas`; the free Wolfram Engine) or give numbers."
                               % (f, x, "" if a is None else " from %s to %s" % (a, b), ", ".join(tried) or "none"))
    finally:
        if extra:
            forget(*extra)
