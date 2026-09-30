import sys; sys.path.insert(0, '.')
from feynsage import *
from feynsage.oneloop import *
# Lee's criterion
kin = Kinematics(['p'], ['s'], {('p','p'): 's'}, euclidean=False)
massless = IntegralFamily('b0', ['l'], kin, [(mom(l=1), 0), (mom(l=1, p=-1), 0)])
massive = IntegralFamily('b1', ['l'], kin, [(mom(l=1), 1), (mom(l=1, p=-1), 1)])
print("massless tadpole sector zero:", massless.is_zero_sector((1,0)), "| massless bubble zero:", massless.is_zero_sector((1,1)),
      "| massive tadpole zero:", massive.is_zero_sector((1,0)))
Q2 = var('Q2')
# the fun phi^4 diagram of the notes: G(1,1) * triangle(1,1,4-D/2)
fun = expand_eps(G(1,1) * triangle_onshell(1, 1, 2 - D/2, 1), 0, loops=2)
print("phi4 two-loop diagram:", fun, "| minus note value:", (fun - (1/(2*eps**2) + 5/(2*eps) + (114 + pi**2)/12)).simplify_full())
tri = expand_eps(triangle_onshell(1,1,1,1), 0)
print("massless triangle, Q2 = 1:", tri)
# FORM: Tr(g_mu g_nu g_rho g_sigma)
t = form.trace(['mu','nu','rho','sigma'])
g = function('g'); mu, nu, rho, sigma = var('mu nu rho sigma')
ref = 4*(g(mu,nu)*g(rho,sigma) - g(mu,rho)*g(nu,sigma) + g(mu,sigma)*g(nu,rho))
print("FORM trace:", t, "| minus textbook:", (t - ref).expand())
