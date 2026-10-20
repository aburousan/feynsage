# Pauli matrices in three dimensions (feynsage.pauli) against explicit 2x2 matrices, and FeynCalc's PauliTrace
# examples (Tr[s^i s^j] = 2 delta^{ij}, Tr[s^i s^j s^k] = 2 i eps^{ijk}, the four-sigma trace) and the identity
# (s.a)(s.b) = a.b + i s.(a x b).
#     sage tests/test_pauli.sage
import itertools, random, time
from feynsage.pauli import *
from feynsage.pauli import DELTA3, EPS3, COMP3, DOT3
t_start = time.time()
FAIL = []
def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAIL.append(name)
sig = [matrix(QQbar, [[0, 1], [1, 0]]), matrix(QQbar, [[0, -I], [I, 0]]), matrix(QQbar, [[1, 0], [0, -1]])]
i, j, k, l = cartesian_indices("i j k l")
a, b, c = vectors3("a b c")
rnd = random.Random(int(3))
vec = {str(v): vector(QQ, [rnd.randint(-9, 9) for _ in range(3)]) for v in (a, b, c)}
def eps(x, y, z):
    return 0 if len({x, y, z}) < 3 else Permutation([x + 1, y + 1, z + 1]).sign()
def num(e, iv):
    """value of an expression in delta3, eps3, comp3, dot3 at index values iv (0, 1, 2); the other (summed)
    indices are summed over"""
    from feynsage.pauli import _monomials, _is_index
    tot = 0
    for m in _monomials(e):
        summed = sorted({str(v) for v in m.variables() if _is_index(v) and str(v) not in iv})
        for vals in itertools.product(range(3), repeat=len(summed)):
            tot += num1(m, dict(iv, **dict(zip(summed, vals))))
    return QQbar(tot)
def num1(e, iv):
    e = SR(e)
    w = [SR.wild(n) for n in range(3)]
    sub = {}
    for F in e.find(DELTA3(w[0], w[1])): x, y = F.operands(); sub[F] = 1 if iv[str(x)] == iv[str(y)] else 0
    for F in e.find(EPS3(w[0], w[1], w[2])): sub[F] = eps(*[iv[str(x)] for x in F.operands()])
    for F in e.find(COMP3(w[0], w[1])): v, x = F.operands(); sub[F] = vec[str(v)][iv[str(x)]]
    for F in e.find(DOT3(w[0], w[1])): u, v = F.operands(); sub[F] = vec[str(u)] * vec[str(v)]
    return QQbar(e.subs(sub))
def explicit(fs, iv):
    M = identity_matrix(QQbar, 2)
    for f in fs:
        M = M * (sig[iv[str(f)]] if str(f) in iv else sum(vec[str(f)][n] * sig[n] for n in range(3)))
    return M
def agrees(fs):
    """symbolic product = explicit product for every value of the free indices; an index written twice is summed
    (in the explicit product as well)"""
    names = [str(f) for f in fs if str(f) in "ijkl"]
    free = sorted({n for n in names if names.count(n) == 1}); summed = sorted({n for n in names if names.count(n) == 2})
    S, V, cc = pauli_simplify(*fs)
    for vals in itertools.product(range(3), repeat=len(free)):
        iv = dict(zip(free, vals))
        mine = num(S, iv) * identity_matrix(QQbar, 2) + sum(num(V, dict(iv, **{str(cc): n})) * sig[n] for n in range(3))
        ex = sum(explicit(fs, dict(iv, **dict(zip(summed, sv)))) for sv in itertools.product(range(3), repeat=len(summed)))
        if mine != ex:
            print("   mismatch for", fs, iv)
            return False
    return True
check("Tr[s^i s^j] = 2 delta^{ij}", bool(pauli_trace(i, j) == 2*DELTA3(i, j)))
check("Tr[s^i s^j s^k] = 2 i eps^{ijk}", bool(pauli_trace(i, j, k) == 2*I*EPS3(i, j, k)))
check("Tr[s^i s^j s^k s^l] = 2 (d^ij d^kl - d^ik d^jl + d^il d^jk)",
      bool((pauli_trace(i, j, k, l) - 2*(DELTA3(i, j)*DELTA3(k, l) - DELTA3(i, k)*DELTA3(j, l) + DELTA3(i, l)*DELTA3(j, k))).expand() == 0))
S, V, cc = pauli_simplify(a, b)
check("(s.a)(s.b) = a.b + i s.(a x b)", bool(S == DOT3(a, b)) and agrees([a, b]))
check("products of up to five factors, indices and vectors mixed, = explicit matrices for every index value",
      all(agrees(fs) for fs in ([i, j, k], [i, i], [a, i, b, j], [i, j, i], [a, b, c, i, a])))
check("Tr[(s.a)(s.b)(s.c)] = 2 i a.(b x c), Tr[s^i s^i] = 6",
      bool(pauli_trace(i, i) == 6) and all(num(pauli_trace(a, b, c), {}) == 2*I*matrix(QQ, [vec['a'], vec['b'], vec['c']]).det() for _ in [0]))
print("time %.0f s" % (time.time() - t_start))
print("ALL PAULI CHECKS PASS" if not FAIL else "FAILED: %s" % FAIL)
