r"""
A scattering process or a decay in one line, written as in a textbook:

    P = process("e- e+ -> mu- mu+")              # Standard Model, tree level, Feynman-'t Hooft gauge
    P.draw()                                     # the diagrams
    P.squared()                                  # |M|^2 averaged over initial spins and colours,
                                                 # summed over final ones, in s, t, u
    P.squared(helicities={"e-": -1, "e+": +1})   # one helicity configuration (CM frame)
    P.dsigma_dcos()                              # d sigma/d cos(theta) in sqrt(s) and cos(theta)
    P.sigma(sqrt_s=91.19, unit="pb")             # a number (Standard-Model values by default)
    process("H -> b b~").width()                 # a decay width

Momenta are p1, p2 (incoming) and p3, p4, ... (outgoing); for 2 -> 2, s = (p1 + p2)^2,
t = (p1 - p3)^2, u = (p1 - p4)^2 and theta is the angle between p1 and p3 in the centre-of-mass
frame.  Particle names: e- e+ mu- mu+ tau- tau+ nu_e nu_e~ ... u u~ d d~ ... t t~ b b~, gamma (or A),
Z, W+, W-, H, g (gluon); an antiparticle is written with ~ (or ubar, dbar, ...).
"""
import itertools
from sage.all import SR, sqrt, pi, I, Integer, factorial, prod

from .qft import tensors as T
from .qft import objects as O
from .qft.operations import spin_sum, polarization_sum
from .qft.color import color_factor
from .topologies import topologies
from .diagrams import insert_fields, conjugate_amplitude
from . import models as MD

GEV2_TO_PB = 0.3893793721e9            # (hbar c)^2 in pb GeV^2

s_, t_, u_ = SR.var('s'), SR.var('t'), SR.var('u')
sqrt_s = SR.var('sqrt_s', latex_name=r'\sqrt{s}', domain='positive')
cos_theta = SR.var('cos_theta', latex_name=r'\cos\theta', domain='real')
_sin = SR.var('fs_sin_theta', latex_name=r'\sin\theta', domain='positive')   # sin(theta) while conjugating
sin_theta = sqrt(1 - cos_theta ** 2)


def _kallen(a, b, c):
    return a ** 2 + b ** 2 + c ** 2 - 2 * a * b - 2 * a * c - 2 * b * c


def _sm_values():
    """Standard-Model numbers (GeV), alpha(0) and on-shell sin(theta_W) = sqrt(1 - MW^2/MZ^2)."""
    MW, MZ = Integer(80379) / 1000, Integer(911876) / 10000
    return {MD.EL: sqrt(4 * pi / Integer(137035999084) * 10 ** 9), MD.MW: MW, MD.MZ: MZ,
            MD.CW: MW / MZ, MD.SW: sqrt(1 - MW ** 2 / MZ ** 2), MD.MH: Integer(12525) / 100,
            MD.GS: sqrt(4 * pi * Integer(118) / 1000),
            MD.ME: Integer(510999) / 10 ** 9, MD.MM: Integer(1056584) / 10 ** 7, MD.ML: Integer(177686) / 10 ** 5,
            MD.MU: Integer(216) / 10 ** 5, MD.MD: Integer(467) / 10 ** 5, MD.MS: Integer(934) / 10 ** 4,
            MD.MC: Integer(127) / 100, MD.MB: Integer(418) / 100, MD.MT: Integer(17276) / 100}


_NAMED = {'alpha': None, 'alpha_s': None, 'sw2': None}


class Process:
    def __init__(self, text, model=None, loops=0, gauge='feynman', exclude_fields=(), massless=(),
                 widths=None):
        if '->' not in text:
            raise ValueError('write the process as "a b -> c d", e.g. "e- e+ -> mu- mu+"')
        lhs, rhs = text.split('->')
        self.text = text.strip()
        self.model = model if model is not None else MD.SM()
        if isinstance(self.model, str):
            self.model = {'SM': MD.SM, 'QED': MD.QED}[self.model]()
        names_in, names_out = lhs.split(), rhs.split()
        if not names_in or not names_out:
            raise ValueError("no particles on one side of %r" % text)
        f = lambda n: self.model.field(n).name
        self.fields_in = [f(n) for n in names_in]
        self.fields_out = [f(n) for n in names_out]
        self.labels = names_in + names_out              # the names as the user wrote them
        self.particles = [self.model.particles[n] for n in self.fields_in + self.fields_out]
        n = len(self.particles)
        self.momenta = list(T.momenta(' '.join('p%d' % i for i in range(1, n + 1))))
        self.p_in, self.p_out = self.momenta[:len(names_in)], self.momenta[len(names_in):]
        if widths:                       # a Goldstone boson propagates with the mass of its vector boson
            widths = dict(widths)
            for v, g in (('Z', 'G0'), ('W+', 'G+'), ('W-', 'G-')):
                if v in widths and g not in widths and gauge != 'unitary':
                    widths[g] = widths[v]
        self.gauge, self.widths, self.loops = gauge, widths, loops
        self._zero = {}
        for name in massless:
            m = self.model.field(name).mass
            if not m.is_trivial_zero():
                self._zero[m] = 0
        self.masses = [p.mass.subs(self._zero) for p in self.particles]
        tops = topologies(loops, len(names_in), len(names_out),
                          exclude=('tadpoles', 'wf') if loops else ())
        self._exclude = exclude_fields
        self.diagrams = insert_fields(tops, self.fields_in, self.fields_out, self.model,
                                      exclude_fields=exclude_fields)
        self.values = _sm_values()
        self.s, self.t, self.u = s_, t_, u_
        self.sqrt_s, self.cos_theta, self.sin_theta = sqrt_s, cos_theta, sin_theta

    # ------------------------------------------------------------------ basics
    def __repr__(self):
        return 'process %s: %d diagram%s' % (self.text, len(self.diagrams), '' if len(self.diagrams) == 1 else 's')

    def draw(self, **kw):
        return self.diagrams.draw(**kw)

    def amplitude(self, diagrams=None):
        """M summed over the diagrams (spinors, polarization vectors and colour left open).
        diagrams: a list of diagram numbers (from 1, as in draw()) to keep only those."""
        return self._amp(diagrams)[0]

    def _amp(self, diagrams=None):
        from .diagrams import _Shared, _sum_amplitudes
        sh = _Shared()
        amps = self.diagrams.amplitudes(self.p_in, self.p_out, gauge=self.gauge, widths=self.widths, shared=sh)
        keep = range(1, len(amps) + 1) if diagrams is None else diagrams
        M = _sum_amplitudes([amps[k - 1] for k in keep if amps[k - 1] is not None])
        return M, list(sh.col.values())

    def _squared_raw(self, diagrams=None):
        M, keep = self._amp(diagrams)
        Mc = conjugate_amplitude(M, keep=keep)
        prod_ = M * Mc
        return prod_ if isinstance(prod_, type(SR(1))) else spin_sum(prod_)

    def _vector_legs(self):
        return [i for i, p in enumerate(self.particles) if p.kind == 'V']

    def _dof(self, p):
        if p.kind == 'F':                # a Standard-Model neutrino has one helicity (no right-handed nu)
            neutrino = p.mass.is_trivial_zero() and SR(p.charge).is_trivial_zero() and not p.colour
            return 1 if neutrino and self.model.name.startswith('SM') else 2
        if p.kind == 'V':
            return 2 if p.mass.is_trivial_zero() else 3
        return 1

    def _colour_dof(self, p):
        return {'3': 3, '3b': 3, '8': 8}.get(p.colour, 1)

    def symmetry_factor(self):
        """1/n! for n identical particles in the final state."""
        out = SR(1)
        for name in set(self.fields_out):
            out /= factorial(self.fields_out.count(name))
        return out

    # ------------------------------------------------------------------ kinematics
    def kinematics(self):
        """The scalar products of the external momenta: p_i^2 = m_i^2, and for 2 -> 2 the dot products
        in s, t, u, for 1 -> 2 in the masses."""
        m2 = [m ** 2 for m in self.masses]
        p = self.momenta
        k = {T.dot(p[i], p[i]): m2[i] for i in range(len(p))}
        if len(self.p_in) == 2 and len(self.p_out) == 2:
            k.update({T.dot(p[0], p[1]): (s_ - m2[0] - m2[1]) / 2, T.dot(p[2], p[3]): (s_ - m2[2] - m2[3]) / 2,
                      T.dot(p[0], p[2]): (m2[0] + m2[2] - t_) / 2, T.dot(p[1], p[3]): (m2[1] + m2[3] - t_) / 2,
                      T.dot(p[0], p[3]): (m2[0] + m2[3] - u_) / 2, T.dot(p[1], p[2]): (m2[1] + m2[2] - u_) / 2})
        elif len(self.p_in) == 1 and len(self.p_out) == 2:
            M2, a, b = m2
            k.update({T.dot(p[1], p[2]): (M2 - a - b) / 2, T.dot(p[0], p[1]): (M2 + a - b) / 2,
                      T.dot(p[0], p[2]): (M2 - a + b) / 2})
        return k

    def _cm(self):
        """Components in the centre-of-mass frame (rest frame for a decay), as functions of sqrt_s and
        cos_theta: p1 along +z, p3 in the x-z plane at angle theta."""
        m2 = [m ** 2 for m in self.masses]
        st, ct = _sin, cos_theta
        if len(self.p_in) == 2 and len(self.p_out) == 2:
            rs = sqrt_s
            E1, E2 = (rs ** 2 + m2[0] - m2[1]) / (2 * rs), (rs ** 2 - m2[0] + m2[1]) / (2 * rs)
            E3, E4 = (rs ** 2 + m2[2] - m2[3]) / (2 * rs), (rs ** 2 - m2[2] + m2[3]) / (2 * rs)
            P = sqrt(_kallen(rs ** 2, m2[0], m2[1])) / (2 * rs)
            Q = sqrt(_kallen(rs ** 2, m2[2], m2[3])) / (2 * rs)
            return [(E1, 0, 0, P), (E2, 0, 0, -P), (E3, Q * st, 0, Q * ct), (E4, -Q * st, 0, -Q * ct)], \
                   [(0, 0), (pi, 0), (None, 0), (None, pi)]
        if len(self.p_in) == 1 and len(self.p_out) == 2:
            M = self.masses[0]
            E2, E3 = (m2[0] + m2[1] - m2[2]) / (2 * M), (m2[0] - m2[1] + m2[2]) / (2 * M)
            Q = sqrt(_kallen(m2[0], m2[1], m2[2])) / (2 * M)
            return [(M, 0, 0, 0), (E2, Q * st, 0, Q * ct), (E3, -Q * st, 0, -Q * ct)], [(0, 0), (None, 0), (None, pi)]
        raise NotImplementedError("helicity amplitudes are made for 2 -> 2 and 1 -> 2")

    # ------------------------------------------------------------------ |M|^2
    def squared(self, average=True, helicities=None, polarization_gauge=None, diagrams=None, nproc=None):
        r"""
        |M|^2.  Without helicities: summed over all spins and polarizations (and colours) and, with
        average=True, divided by the number of initial spin and colour states (2 for a massive or
        charged fermion, 1 for a neutrino, which is only left-handed); the result is in
        s, t, u for 2 -> 2 (in the masses for 1 -> 2).  Photons are summed with -g, gluons with the
        physical sum (reference vector: another external momentum), massive vectors with -g + kk/M^2.

        helicities={"e-": -1, "e+": +1, "W+": 0, "W-": +1}: one helicity configuration in the
        centre-of-mass frame (rest frame for a decay), as a function of sqrt_s and cos_theta.
        Fermions: +1 or -1 (twice the helicity; "R"/"L" also accepted), only for massless fermions
        (make them massless with process(..., massless=["e"])).  Vectors: +1, -1 or 0.
        Particles not named are summed over; nothing is averaged.
        diagrams=[1, 3]: only these diagrams (numbers as in draw()).
        nproc: worker processes (feynsage.parallel; None = all cores, 1 = serial).  Without helicities
        the products M_i M_j* of pairs of diagrams are worked out in parallel.
        A particle that appears twice is named by its position: {3: +1} is the third particle.
        """
        if helicities:
            return self._helicity(helicities, diagrams)
        r = self._unpolarized(diagrams, polarization_gauge, nproc)
        if average:
            r = r / prod(self._dof(p) * self._colour_dof(p) for p in self.particles[:len(self.p_in)])
        return r

    def _unpolarized(self, diagrams, gauge, nproc):
        from .parallel import pmap, resolve
        from .diagrams import _Shared
        sh = _Shared()
        amps = self.diagrams.amplitudes(self.p_in, self.p_out, gauge=self.gauge, widths=self.widths, shared=sh)
        keep = list(sh.col.values())
        idx = [k for k in (diagrams or range(1, len(amps) + 1)) if amps[k - 1] is not None]
        kin, legs = self.kinematics(), self._vector_legs()

        def finish(prod_):
            r = prod_ if isinstance(prod_, type(SR(1))) else spin_sum(prod_)
            r = self._colour(r)                  # colour first: it merges many terms
            r = self._sum_vectors(r, legs, gauge, kin)
            return self._no_eps(SR(r).subs(kin)).subs(self._zero)
        if resolve(nproc) > 1 and len(idx) >= 2:
            conj = {k: conjugate_amplitude(amps[k - 1], keep=keep) for k in idx}
            pairs = [(i, j) for i in idx for j in idx]
            return sum(pmap(lambda ij: finish(amps[ij[0] - 1] * conj[ij[1]]), pairs, nproc=nproc), SR(0))
        from .diagrams import _sum_amplitudes
        M = _sum_amplitudes([amps[k - 1] for k in idx])
        return finish(M * conjugate_amplitude(M, keep=keep))

    def _sum_vectors(self, r, legs, gauge=None, kin=None, dim=4):
        """Polarization sums over these legs; kin (the scalar products) is put in after every sum,
        which keeps the expression small (four gluons: minutes -> seconds)."""
        if kin:
            r = SR(r).subs(kin)
        for i in legs:
            p = self.particles[i]
            e = SR.var('eps%d' % (i + 1))
            if str(e) not in T._POL:
                continue
            if not p.mass.is_trivial_zero():
                r = polarization_sum(r, e, dim=dim)
            elif p.colour == '8' or gauge == 'physical':
                ref = [q for j, q in enumerate(self.momenta) if j != i][0]
                ref = next((q for j, q in enumerate(self.momenta) if j != i and self.particles[j].mass.is_trivial_zero()), ref)
                r = polarization_sum(r, e, gauge='physical', n=ref, dim=dim)
            else:
                r = polarization_sum(r, e, gauge='feynman', dim=dim)
            if kin:
                r = SR(r).subs(kin)
        return r

    def _no_eps(self, r):
        """With at most four external momenta only three are independent, so every epsilon tensor of
        external momenta alone is zero (p4 = p1 + p2 - p3)."""
        if len(self.momenta) > 4:
            return r
        w = [SR.wild(k) for k in range(4)]
        mom = set(self.momenta)
        hits = [x for x in r.find(T.EPS(*w)) if all(a in mom for a in x.operands())]
        return r.subs({x: 0 for x in hits}) if hits else r

    def _colour(self, r):
        if any(p.colour for p in self.particles) or any(self.model.particles[f].colour for d in self.diagrams
                                                       for f, _, _ in d.fields):
            r = color_factor(r, N=3)
        return r

    def _leg(self, key):
        if not isinstance(key, str):
            return int(key) - 1
        hits = [i for i, l in enumerate(self.labels) if l == key]
        if not hits:
            name = self.model.field(key).name
            hits = [i for i, p in enumerate(self.particles) if p.name == name]
        if len(hits) != 1:
            raise ValueError("%r is not one particle of %s; name it by its position (1, 2, ...)" % (key, self.text))
        return hits[0]

    def _helicity(self, hels, diagrams=None):
        key, epsc, allc = self._hel_parse(hels, diagrams)
        cache = self.__dict__.setdefault('_hcache', {})
        if key not in cache:
            cache[key] = self._hel_trace(key, allc)
        return _tidy(_frame(cache[key], {**allc, **epsc}).subs(self._zero))

    def helicity_table(self, configs, nproc=None):
        """|M|^2 for many helicity configurations (a list of dicts as for squared(helicities=...)).
        The different fermion-helicity traces are made in parallel; returns a list of results."""
        from .parallel import pmap
        parsed = [self._hel_parse(h, None) for h in configs]
        cache = self.__dict__.setdefault('_hcache', {})
        todo = sorted({k for k, _, _ in parsed if k not in cache}, key=str)
        if todo:
            allc = parsed[0][2]
            for k, r in zip(todo, pmap(lambda k: self._hel_trace(k, allc), todo, nproc=nproc, probe=False)):
                cache[k] = r
        return pmap(lambda kv: _tidy(_frame(cache[kv[0]], {**kv[2], **kv[1]}).subs(self._zero)), parsed, nproc=nproc)

    def _hel_parse(self, hels, diagrams):
        comps, angles = self._cm()
        basis = T.momenta('fs_e0 fs_ex fs_ey fs_ez')
        bcomp = {basis[0]: (1, 0, 0, 0), basis[1]: (0, 1, 0, 0), basis[2]: (0, 0, 1, 0), basis[3]: (0, 0, 0, 1)}
        allc = dict(bcomp)
        allc.update({self.momenta[i]: comps[i] for i in range(len(comps))})
        vec = lambda c: sum(SR(x) * b for x, b in zip(c, basis))
        hands, vals, labelled = [], {}, []
        for key, h in hels.items():
            i = self._leg(key)
            p = self.particles[i]
            if isinstance(h, str):
                h = {'R': 1, 'L': -1, '+': 1, '-': -1, '0': 0}[h.upper()]
            if p.kind == 'F':
                if not p.mass.subs(self._zero).is_trivial_zero():
                    raise ValueError("helicity of a massive fermion (%s): make it massless with "
                                     "process(..., massless=[%r])" % (self.labels[i], self.labels[i]))
                particle = not p.name.endswith('~')
                hands.append((i, 'R' if (h > 0) == particle else 'L'))
            elif p.kind == 'V':
                E, kx, ky, kz = comps[i]
                kk = sqrt(kx ** 2 + ky ** 2 + kz ** 2)
                th, ph = angles[i]
                if th is None:
                    cth, sth = (cos_theta, _sin) if ph == 0 else (-cos_theta, _sin)
                else:
                    cth, sth = SR(th).cos(), SR(th).sin()
                cph, sph = SR(ph).cos(), SR(ph).sin()
                if h == 0:
                    if p.mass.is_trivial_zero():
                        raise ValueError("a massless vector has no helicity 0")
                    eps = (kk / p.mass, E * sth * cph / p.mass, E * sth * sph / p.mass, E * cth / p.mass)
                else:
                    tht = (0, cth * cph, cth * sph, -sth)
                    phh = (0, -sph, cph, 0)
                    eps = tuple(-h * (a + h * I * b) / sqrt(2) for a, b in zip(tht, phh))
                # eps_i is the vector eps, eps_i_c its conjugate (M has eps for an incoming vector, eps* for an outgoing one)
                vals[SR.var('eps%d' % (i + 1))] = tuple(eps)
                vals[SR.var('eps%d_c' % (i + 1))] = tuple(x.conjugate() for x in eps)
                labelled.append(i)
            else:
                raise ValueError("%s has no helicity" % self.labels[i])
        key = (tuple(sorted(hands)), tuple(sorted(labelled)), tuple(diagrams) if diagrams else None)
        return key, vals, allc

    def _hel_trace(self, key, allc):
        """One trace for this choice of fermion helicities, the labelled polarization vectors kept as
        vectors (their components are put in by _frame)."""
        hands, labelled, diagrams = key
        M, keep = self._amp(list(diagrams) if diagrams else None)
        for i, hand in hands:
            M = O.chiral(M, self.momenta[i], hand)
        # the polarization vectors stay vectors here; their components go in afterwards (_frame)
        prod_ = M * conjugate_amplitude(M, keep=keep)
        r = prod_ if isinstance(prod_, type(SR(1))) else SR(spin_sum(prod_))
        r = self._colour(r)
        r = self._sum_vectors(r, [i for i in self._vector_legs() if i not in labelled])
        return _frame(r, allc)

    # ------------------------------------------------------------------ one loop
    def loop_diagrams(self):
        """The one-loop diagrams of the process (no tadpoles, no corrections on external legs)."""
        if not hasattr(self, '_loop_diags'):
            n_in, n_out = len(self.fields_in), len(self.fields_out)
            self._loop_diags = insert_fields(topologies(1, n_in, n_out, exclude=('tadpoles', 'wf')),
                                             self.fields_in, self.fields_out, self.model,
                                             exclude_fields=self._exclude)
        return self._loop_diags

    def virtual(self, average=True, xi=1, nproc=None, diagrams=None):
        r"""
        The one-loop correction to |M|^2: sum over spins, polarizations and colours of M0^* M1, with
        M0 the tree amplitude and M1 the unrenormalized one-loop amplitude (all diagrams of
        loop_diagrams(), Feynman-'t Hooft gauge for xi = 1), averaged like squared().  Everything is in
        d = 4 - 2 eps dimensions (Dirac algebra, polarization sums); the result is a Laurent series in
        eps with A0, B0, C0, D0 written out (finite_part, pole_parts).  gamma5 anticommutes (NDR); with at
        most four external legs the traces with one gamma5 are dropped, as their epsilon tensors vanish
        after the integration.  The virtual correction to the
        cross section is 2 Re of it; for 2 -> 2 the result is in s and t (u = sum m^2 - s - t).
        diagrams=[1, 4]: only those loop diagrams.  The diagrams are done in parallel (nproc).
        """
        from .diagrams import _Shared, _loop_integrand, _to_loop_string, _rational_poly
        from .qft.operations import contract
        from .pv import loop
        from .parallel import pmap
        sh = _Shared()
        amps = self.diagrams.amplitudes(self.p_in, self.p_out, gauge='feynman', shared=sh)
        from .diagrams import _sum_amplitudes
        M0 = _sum_amplitudes([a for a in amps if a is not None])
        keep = list(sh.col.values())
        M0c = conjugate_amplitude(M0, keep=keep)
        lmom = T._declare('l', 'vector')
        lmom = lmom[0] if isinstance(lmom, (list, tuple)) else lmom
        kin = self.kinematics()
        if len(self.p_in) == 2 and len(self.p_out) == 2:     # u = sum m^2 - s - t, so the invariants are on shell
            usub = {u_: sum(m ** 2 for m in self.masses) - s_ - t_}
            kin = {k_: SR(v_).subs(usub) for k_, v_ in kin.items()}
        kin_str = {}
        for key, val in kin.items():
            a, b = key.operands()
            kin_str[('%s^2' % a) if a == b else ('%s.%s' % (a, b))] = str(SR(val).subs(self._zero))
        legs = self._vector_legs()
        loops = self.loop_diagrams()
        idx = diagrams or range(1, len(loops) + 1)
        jobs = []
        for k in idx:
            terms, ext, col = _loop_integrand(loops[k - 1], self.p_in, self.p_out, lmom, xi=SR(xi),
                                              prefix='v%d' % k, shared=sh, spinors=True, conserve=True)
            if isinstance(col, tuple):
                raise NotImplementedError("four-gluon vertices in virtual(); use loop_amplitude")
            jobs.append((terms, ext, col))
        n_in = len(self.p_in)
        pol = {}
        for i in legs:
            e = SR.var('eps%d' % (i + 1))
            pol[i] = e if i < n_in else T.conjugate_vector(e)

        def one(job):
            terms, ext, col = job
            total = SR(0)
            for num, props in terms:
                x = col * num
                for i in legs:                       # the external vectors: their polarization vectors
                    x = x * T.comp(pol[i], ext[i])
                # with at most four external momenta (three independent) every epsilon tensor vanishes after
                # the loop integration, so the traces with one gamma5 are dropped (anticommuting gamma5)
                g5 = 'NDR-even' if len(self.momenta) <= 4 else 'NDR'
                x = spin_sum(x * M0c, dim='d', gamma5_scheme=g5) if not isinstance(x * M0c, type(SR(1))) else x * M0c
                x = self._colour(x)
                x = self._sum_vectors(x, legs, None, kin, dim='d')
                x = SR(contract(SR(x).subs(kin), dim='d')).subs(kin).subs(self._zero).expand()
                if not props:
                    total += -I * x
                    continue
                pr, back = [], {}
                for q, m in props:
                    qs = SR(q)
                    if qs.coefficient(lmom) == -1:
                        qs = -qs
                    m = SR(m).subs(self._zero)
                    if not _rational_poly(m):
                        if m not in back.values():
                            back[SR.var('fsmass%d' % len(back))] = m
                        m = [kk for kk, vv in back.items() if (vv - m).is_trivial_zero()][0]
                    pr.append([str(qs), str(m)])
                res = loop(_to_loop_string(x, 'l'), *pr, kin=kin_str)
                val = res.coefficients().get('1', SR(0)) if hasattr(res, 'coefficients') else SR(res)
                if back:
                    val = SR(val).subs(back)
                total += val / (16 * pi ** 2)
            return total
        r = sum(pmap(one, jobs, nproc=nproc), SR(0))
        if average:
            r = r / prod(self._dof(p) * self._colour_dof(p) for p in self.particles[:len(self.p_in)])
        return r

    # ------------------------------------------------------------------ rates
    def width(self):
        """Gamma for a 1 -> 2 decay: |p| /(8 pi M^2) |M|^2 (averaged), times 1/2 for identical particles."""
        if len(self.p_in) != 1 or len(self.p_out) != 2:
            raise NotImplementedError("width() is for 1 -> 2 decays")
        M2, a, b = [m ** 2 for m in self.masses]
        P = sqrt(_kallen(M2, a, b)) / (2 * self.masses[0])
        return self.symmetry_factor() * P / (8 * pi * M2) * self.squared()

    def dsigma_dcos(self, helicities=None):
        """d sigma/d cos(theta) for 2 -> 2 in sqrt_s and cos_theta (GeV^-2):
        |M|^2 /(32 pi s) |p_f|/|p_i|, times 1/n! for identical final particles."""
        if len(self.p_in) != 2 or len(self.p_out) != 2:
            raise NotImplementedError("dsigma_dcos() is for 2 -> 2")
        m2 = [m ** 2 for m in self.masses]
        s = sqrt_s ** 2
        Pi = sqrt(_kallen(s, m2[0], m2[1])) / (2 * sqrt_s)
        Pf = sqrt(_kallen(s, m2[2], m2[3])) / (2 * sqrt_s)
        if helicities:
            M2 = self.squared(helicities=helicities)
        else:
            M2 = self.squared().subs(self.to_angles())
        return self.symmetry_factor() * M2 / (32 * pi * s) * Pf / Pi

    def to_angles(self):
        """s, t, u in sqrt_s and cos_theta (2 -> 2, CM frame)."""
        comps, _ = self._cm()
        mdot = lambda a, b: a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
        d = lambda a, b: tuple(x - y for x, y in zip(a, b))
        return {s_: sqrt_s ** 2, t_: _tidy(mdot(d(comps[0], comps[2]), d(comps[0], comps[2]))).simplify_full(),
                u_: _tidy(mdot(d(comps[0], comps[3]), d(comps[0], comps[3]))).simplify_full()}

    def numeric(self, expr, **values):
        """expr with the Standard-Model numbers put in.  Keywords: sqrt_s=, cos_theta=, and any
        parameter by name (MW=80.4, MZ=..., MH=..., ME=..., alpha=1/128, alpha_s=0.12, sw2=0.23)."""
        v = dict(self.values)
        sub = {}
        for k, x in values.items():
            if k == 'alpha':
                v[MD.EL] = sqrt(4 * pi * SR(x))
            elif k == 'alpha_s':
                v[MD.GS] = sqrt(4 * pi * SR(x))
            elif k == 'sw2':
                v[MD.SW], v[MD.CW] = sqrt(SR(x)), sqrt(1 - SR(x))
            elif k == 'sqrt_s':
                sub[sqrt_s] = x
            elif k == 'cos_theta':
                sub[cos_theta] = x
            else:
                sym = getattr(MD, k, None)
                if sym is None:
                    raise KeyError("unknown parameter %r" % k)
                v[sym] = x
        if 'MW' in values and 'sw2' not in values:          # keep M_W = M_Z c_W
            v[MD.CW] = v[MD.MW] / v[MD.MZ]
            v[MD.SW] = sqrt(1 - v[MD.CW] ** 2)
        if 'MZ' in values and 'sw2' not in values:
            v[MD.CW] = v[MD.MW] / v[MD.MZ]
            v[MD.SW] = sqrt(1 - v[MD.CW] ** 2)
        return SR(expr).subs(sub).subs(v)

    def sigma(self, sqrt_s, unit='GeV^-2', digits=15, helicities=None, **values):
        """The total cross section at this sqrt(s) (numerical integral over cos(theta))."""
        import mpmath
        from sage.all import fast_callable, CDF
        f = self.numeric(self.dsigma_dcos(helicities=helicities), sqrt_s=sqrt_s, **values)
        g = fast_callable(SR(f), vars=[cos_theta], domain=CDF)
        val = float(mpmath.quad(lambda c: float(g(float(c)).real()), [-1, 0, 1]))
        return val * GEV2_TO_PB if unit == 'pb' else val


def _tidy(expr):
    """sin(theta) back as sqrt(1 - cos^2), sin^2 as 1 - cos^2."""
    e = SR(expr).subs({_sin ** 2: 1 - cos_theta ** 2})
    return e.subs({_sin: sqrt(1 - cos_theta ** 2)})


def _frame(expr, comps):
    """Every Dot and Epsilon of vectors with known components replaced by its value."""
    mdot = lambda a, b: a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    from sage.all import matrix
    w = [SR.wild(k) for k in range(4)]
    e = SR(expr)
    sub = {}
    for x in e.find(T.DOT(w[0], w[1])):
        a, b = x.operands()
        if a in comps and b in comps:
            sub[x] = mdot(comps[a], comps[b])
    for x in e.find(T.EPS(*w)):
        args = x.operands()
        if all(a in comps for a in args):
            sub[x] = matrix(SR, [[comps[a][0], -comps[a][1], -comps[a][2], -comps[a][3]] for a in args]).det()
    return e.subs(sub)


def process(text, model=None, **kw):
    return Process(text, model=model, **kw)


process.__doc__ = __doc__ + """
Methods of the result P: P.draw(), P.amplitude(), P.squared(), P.squared(helicities={...}),
P.helicity_table(...), P.dsigma_dcos(), P.sigma(sqrt_s), P.width(), P.numeric(expr, ...), P.kinematics(),
P.to_angles(), P.loop_diagrams(), P.virtual() (the one-loop correction sum M1 M0^*)."""
