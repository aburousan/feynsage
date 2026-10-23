# integrate_checked: several backends, every answer checked numerically; and the places that use it
# (integrate_sectors(exact=True), principal_value, integrate_termwise).
#     sage tests/test_integrators.sage
# With a Wolfram kernel or FriCAS installed those backends are tested too; without them they are skipped.
import time
from feynsage import *
from feynsage.integrators import _has
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
def raises(f, err):
    try:
        f()
    except err:
        return True
    return False
def zero(e, digits=30):
    return abs(CC(SR(e).n(digits=digits + 10))) < 10^(-digits)

x, a, s, MZ = var('x a s MZ')
print("integrators here:", available_integrators())
cases = [("Int_0^1 log(1 + x)/x = pi^2/12", log(1 + x)/x, 0, 1, (), pi^2/12),
         ("Int_0^1 x log(a + x), a > 0", x*log(a + x), 0, 1, (a > 0,),
          a/2 - 1/4 + a^2/2*log(a) - (a^2 - 1)/2*log(1 + a)),
         ("Int_0^oo dx/(x^2 + a^2) = pi/(2a), a > 0", 1/(x^2 + a^2), 0, oo, (a > 0,), pi/(2*a)),
         ("Int_0^1 log(x) log(1 - x) = 2 - pi^2/6", log(x)*log(1 - x), 0, 1, (), 2 - pi^2/6)]
backends = [b for b in ("maxima", "sympy", "giac", "fricas", "wolfram") if _has(b)]
for name, f, lo, hi, asm, ref in cases:
    for b in [None] + backends:
        try:
            r, who = integrate_checked(f, x, lo, hi, assume_=asm, info=True, integrators=None if b is None else [b])
        except IntegrationError:
            if b is not None:
                print("     (%s alone: no checked answer for %s)" % (b, name))
                continue
            check(name + " (default order)", False)
            continue
        pt = {a: 7/3} if asm else {}
        check("%s with %s" % (name, who if b is None else b), zero((r - ref).subs(pt)))

# an antiderivative is checked by differentiating it
check("Int x exp(-a x^2) dx = -exp(-a x^2)/(2a) (indefinite, checked by d/dx)",
      bool((integrate_checked(x*exp(-a*x^2), x) + exp(-a*x^2)/(2*a)).simplify_full() == 0))

# polynomials need no backend: d sigma/d cos theta of e+e- -> mu+mu- with the Z width.  Maxima alone, with M_Z
# still a symbol in the complex Breit-Wigner factors, returns a wrong sigma (14.7 pb instead of 1754.6 pb at the peak)
Pw = process("e- e+ -> mu- mu+", massless=["e", "mu"], widths={"Z": 2.4952})
c = Pw.cos_theta
vals = dict(alpha=1/137.036, MZ=91.1876, sw2=0.23122)
r, who = integrate_checked(Pw.dsigma_dcos(), c, -1, 1, info=True)
from feynsage.process import GEV2_TO_PB
check("sigma(e+e- -> mu+mu-) at the Z peak by the polynomial rule = Pw.sigma (1754.5568 pb)",
      who == "polynomial" and abs(GEV2_TO_PB*CC(Pw.numeric(r, sqrt_s=91.1876, **vals)).real()
                                  - Pw.sigma(91.1876, unit="pb", **vals)) < 1e-6)
check("Maxima alone on the same integral does not pass the check",
      raises(lambda: integrate_checked(Pw.dsigma_dcos(), c, -1, 1, integrators=["maxima"]), IntegrationError))

# Giac sometimes returns the integral unevaluated; that is never accepted as an answer
check("an unevaluated integral is not an answer: Giac on Int_0^1 log(1 + x)/x",
      raises(lambda: integrate_checked(log(1 + x)/x, x, 0, 1, integrators=["giac"]), IntegrationError))

# principal value and term-by-term integration
check("PV Int_0^3 x/(x - 1) dx = 3 + log 2", zero(principal_value(x/(x - 1), x, 0, 3, 1) - 3 - log(2)))
assume(a > 0)
pv2 = principal_value(1/((x - 1)*(x + a)), x, 0, 2, 1)
forget(a > 0)
check("PV Int_0^2 dx/((x - 1)(x + a)) = -log((a + 2)/a)/(a + 1) for symbolic a > 0 (checked inside at random a)",
      zero((pv2 + log((a + 2)/a)/(a + 1)).subs(a=7/3)) and zero((pv2 + log((a + 2)/a)/(a + 1)).subs(a=40)))
# two poles inside: both must be given (this once passed with only c = b, and its check integrated over the other pole)
bp = var('bp'); assume(bp > 0, bp < 1)
check("PV Int_-1^1 dx/(b^2 - x^2) with poles [-b, b] = 2 arctanh(b)/b",
      zero((principal_value(1/(bp^2 - x^2), x, -1, 1, [-bp, bp]) - 2*arctanh(bp)/bp).subs(bp=1/3)))
check("a pole left out of c is an error, not a wrong number",
      raises(lambda: principal_value(1/(bp^2 - x^2), x, -1, 1, bp), IntegrationError))
forget(bp > 0, bp < 1)
check("integrate_termwise(x + x log x, x, 0, 1) = 1/4", zero(integrate_termwise(x + x*log(x), x, 0, 1) - 1/4))

# exact sector decomposition (this used to stop with Maxima's "limit: direction must be either 'plus' or 'minus'")
x1, x2, x3, x4, t = var('x1 x2 x3 x4 t')
secs = sector_decompose([(x1 + x2 + x3, -1 + 2*eps), (-s*x1*x3, -1 - eps)], [x1, x2, x3])
tri = (laurent(-gamma(1 + eps)*exp(eps*euler_gamma), eps, 2)*integrate_sectors(secs, eps, 0, {s: -1}, exact=True)).expand()
exact = laurent(-exp(eps*euler_gamma)*gamma(1 + eps)*gamma(-eps)^2/gamma(1 - 2*eps), eps, 0)
check("massless triangle, exact sector integration = Gamma-function formula (40 digits)",
      all(zero((tri - exact).expand().coefficient(eps, k), 40) for k in (-2, -1, 0)))
secsB = sector_decompose([(x1 + x2 + x3 + x4, 2*eps), (-s*x1*x3 - t*x2*x4, -2 - eps)], [x1, x2, x3, x4])
box = (laurent(gamma(2 + eps)*exp(eps*euler_gamma), eps, 2)*integrate_sectors(secsB, eps, 0, {s: -1, t: -1}, exact=True)).expand()
check("massless box at s = t = -1, exact sector integration = 4/eps^2 - 4 pi^2/3 (40 digits)",
      all(zero((box - (4/eps^2 - 4*pi^2/3)).expand().coefficient(eps, k), 40) for k in (-2, -1, 0)))

print("time %.0f s" % (time.time() - t_start))
print("ALL INTEGRATOR CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
