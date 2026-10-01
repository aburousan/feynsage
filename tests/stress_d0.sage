# Long randomized check (run on a server): random real kinematics, zero to four massless lines.
# For each IR-finite point: d0_value under 6 random relabellings of the lines must agree, and for
# every 4th point d0_closed (exact expression, checked against the number internally) is built.
import sys, time, random, itertools; sys.path.insert(0, '.')
from feynsage.scalar import d0_value, d0_closed
from feynsage.ir import classify_d0
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 300
def relabel(a, pm):
    P = {(0, 1): a[0], (1, 2): a[1], (2, 3): a[2], (0, 3): a[3], (0, 2): a[4], (1, 3): a[5]}
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    return [pp(pm[0], pm[1]), pp(pm[1], pm[2]), pp(pm[2], pm[3]), pp(pm[0], pm[3]), pp(pm[0], pm[2]), pp(pm[1], pm[3])] + [a[6 + q] for q in pm]
t0, bad, done, skipped = time.time(), 0, 0, 0
perms = list(itertools.permutations(range(4)))
while done < npts:
    sc = random.choice([1, 4, 16])
    inv = [QQ(random.randint(-8 * sc, 8 * sc)) / 8 for _ in range(6)]
    ms = [QQ(random.randint(1, 12)) / 4 for _ in range(4)]
    for j in random.sample(range(4), random.choice([0, 0, 1, 2, 3, 4])):
        ms[j] = 0
    a = inv + ms
    try:
        if classify_d0(*a) is not None:
            continue
    except NotImplementedError:
        continue
    try:
        v0 = d0_value(*a)
        vals = [d0_value(*relabel(a, pm)) for pm in random.sample(perms, 6)]
        if done % 4 == 0:
            d0_closed(*a)
    except ZeroDivisionError:
        skipped += 1; done += 1; continue
    except Exception as e:
        bad += 1; done += 1; print("ERROR", a, repr(e)[:120]); continue
    sp = max(abs(v - v0) for v in vals) / max(1e-300, abs(v0))
    if sp > 1e-15:
        bad += 1; print("INCONSISTENT %.1e" % sp, a, v0, vals[:2])
    done += 1
print("stress: %d points, %d bad, %d exceptional (no regular labelling), %.0f s" % (done, bad, skipped, time.time() - t0))
