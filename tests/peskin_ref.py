"""Peskin and Schroeder (7.91) evaluated with mpmath, as an independent reference.
Plain Python on purpose: the Sage preparser would turn the literals into Sage numbers,
which older mpmath versions reject."""
import mpmath


def pi_hat(sv, mv=1.0):
    """Pi_hat(q^2 = sv) in units of alpha/(4 pi): -8 Int x(1-x) log(m^2/(m^2 - x(1-x) q^2 - i0))."""
    sv, mv = mpmath.mpf(float(sv)), mpmath.mpf(float(mv))
    tiny = mpmath.mpf('1e-40')
    f = lambda x: x * (1 - x) * (mpmath.log(mv ** 2) - mpmath.log(mpmath.mpc(mv ** 2 - x * (1 - x) * sv, -tiny)))
    pts = [0, 0.5, 1]
    if sv > 4 * mv ** 2:                     # split at the zeros of m^2 - x(1-x) s (log singularities)
        r = mpmath.sqrt(1 - 4 * mv ** 2 / sv)
        pts = [0, (1 - r) / 2, 0.5, (1 + r) / 2, 1]
    with mpmath.workdps(30):
        return complex(-8 * mpmath.quad(f, pts))
