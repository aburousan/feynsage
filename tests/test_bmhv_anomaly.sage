# The Adler-Bell-Jackiw anomaly from the triangle diagrams in dimensional regularisation with the BMHV gamma5,
# as in Peskin and Schroeder, Sec. 19.2, eqs. (19.47)-(19.60): i q_mu M^{mu nu lambda} (both diagrams) must be
# (e^2/(2 pi^2)) eps^{alpha lambda beta nu} k_alpha p_beta.  It comes only from the l_perp^2 (hat) terms: (d - 4)
# times a 1/(d - 4) pole, eqs. (19.55)-(19.58).  Off-shell photons (p^2, k^2 symbols), so there are no IR poles.
#     sage tests/test_bmhv_anomaly.sage
import time
from feynsage import *
from feynsage.bmhv import bmhv_trace, bmhv_loop, canon_epsilon
from feynsage.pv import eps
t_start = time.time()
p, k, l = momenta("p k l"); lam, nu = lorentz_indices("lam nu")
kin = {"p^2": "P2", "k^2": "K2", "p.k": "PK"}
q = p + k
# (19.49): i q_mu x diagram 1 = e^2 Int d^dl/(2 pi)^d tr[q/ g5 (l/ - k/) g^lam l/ g^nu (l/ + p/)] / ((l - k)^2 l^2 (l + p)^2);
# diagram 2: (p, nu) <-> (k, lam).  Int d^dl/(2 pi)^d = i/(16 pi^2) x loop() (+ O(eps), harmless: the total is finite)
total = SR(0)
for (a, A), (b, B) in [((k, lam), (p, nu)), ((p, nu), (k, lam))]:
    N = bmhv_trace(slash(q)*gamma5()*slash(l - a)*gamma(A)*slash(l)*gamma(B)*slash(l + b), four_dim=(p, k))
    total += bmhv_loop(N, l, [["l - %s" % a, 0], ["l", 0], ["l + %s" % b, 0]], kin, four_dim=(p, k), four_indices=(lam, nu))
e_ = var('e_')
res = (e_^2 * I/(16*pi^2) * total).expand()
poles = [res.coefficient(eps, -n) for n in (1, 2)]
fin = res.coefficient(eps, 0).subs({dot(p, p): var('P2'), dot(k, k): var('K2'), dot(p, k): var('PK')})
ref = canon_epsilon(e_^2/(2*pi^2) * epsilon(k, lam, p, nu))
poles = [SR(x).simplify_rational() for x in poles]
diff = (fin - ref).simplify_rational()
print("1/eps and 1/eps^2 parts:", poles)
print("finite part:", canon_epsilon(fin.simplify_rational()), "  minus (e^2/(2 pi^2)) eps^{alpha lambda beta nu} k_alpha p_beta:", diff)
ok = all(x.is_zero() for x in poles) and diff.is_zero()
print("time %.0f s" % (time.time() - t_start))
print("BMHV ANOMALY CHECK PASSES" if ok else "FAILED")
