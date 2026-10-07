r"""
Models for diagram generation: particles, propagators and vertices.

    sm = SM()                 # the Standard Model with QCD, as FeynArts' {"SM", "SMQCD"}
    qed = QED()               # electrons, muons, taus and the photon

The Standard Model is the tree-level part of FeynArts' SM.mod and SMQCD.mod (H. Eck, S. Kueblbeck,
T. Hahn), which follow A. Denner, Fortschr. Phys. 41 (1993) 307: the same fields, the same signs of
the couplings, R_xi gauge with xi = 1 (Feynman-'t Hooft) by default and no quark mixing.

Fields are named by strings.  A particle and its antiparticle: 'e' and 'e~', 'W-' and 'W+', 'G-' and
'G+'; self-conjugate fields have one name ('A', 'Z', 'H', 'G0', 'g').  Every vertex lists its fields
as all incoming, in the order of the model file, with a Lorentz structure and coefficients:

    FFV  (psibar, psi, V):  gamma^mu3 (cL P_L + cR P_R)
    FFS  (psibar, psi, S):  cL P_L + cR P_R
    VVV:  g^{12} (k2 - k1)^3 + g^{23} (k3 - k2)^1 + g^{31} (k1 - k3)^2
    VVVV: c1 g^{12} g^{34} + c2 g^{13} g^{24} + c3 g^{14} g^{32}
    SSV:  (k1 - k2)^3        SVV: g^{23}        SSVV: g^{34}        SSS, SSSS, SUU: 1
    UUV:  c1 k1^3 + c2 k2^3
all momenta incoming.  The vertex factor (Feynman rule, with its i) is factor * structure.
"""
from sage.all import SR, sqrt, I


class Particle:
    def __init__(self, name, anti, kind, mass=0, charge=0, colour=None, latex=None, xi=None):
        self.name, self.anti, self.kind = name, anti, kind          # kind: F, V, S, U
        self.mass, self.charge, self.colour = SR(mass), charge, colour
        self.latex = latex or name
        self.xi = xi                                                # name of the gauge parameter

    @property
    def self_conjugate(self):
        return self.name == self.anti

    def __repr__(self):
        return self.name


class Vertex:
    def __init__(self, fields, kind, factor, coeffs=(1,), colour=None):
        self.fields, self.kind = tuple(fields), kind
        self.factor, self.coeffs = SR(factor), tuple(SR(c) for c in coeffs)
        self.colour = colour            # None, 'delta' (q qbar), 'T' (q qbar g), 'f' (ggg), 'ffff' (gggg)

    def __repr__(self):
        return 'Vertex(%s)' % ', '.join(self.fields)


class Model:
    def __init__(self, name, particles, vertices, params=None):
        self.name = name
        self.particles = {p.name: p for p in particles}
        for p in particles:                                         # antiparticle entries
            if p.anti not in self.particles:
                self.particles[p.anti] = Particle(p.anti, p.name, p.kind, p.mass, -p.charge,
                                                  _anti_colour(p.colour), _anti_latex(p.latex), p.xi)
        self.vertices = [v for v in vertices if any(not SR(c).is_trivial_zero() for c in v.coeffs)
                         and not SR(v.factor).is_trivial_zero()]
        self.params = params or {}
        self.aliases = {}
        self._by_set = {}
        for v in self.vertices:
            self._by_set.setdefault(tuple(sorted(v.fields)), []).append(v)

    def field(self, name):
        name = self.aliases.get(name, name)
        if name not in self.particles:
            raise KeyError("no field %r in model %s (fields: %s)" % (name, self.name, ', '.join(sorted(self.particles))))
        return self.particles[name]

    def anti(self, name):
        return self.field(name).anti

    def find(self, fields):
        """The vertices whose field content is this multiset (list of names, all incoming)."""
        return self._by_set.get(tuple(sorted(fields)), [])

    def __repr__(self):
        return 'Model %s: %d fields, %d vertices' % (self.name, len(self.particles), len(self.vertices))


def _anti_latex(t):
    for a, b in (('^-', '^+'), ('_-', '_+')):
        if a in t:
            return t.replace(a, b)
        if b in t:
            return t.replace(b, a)
    return r'\bar{%s}' % t


def _anti_colour(c):
    return {'3': '3b', '3b': '3', '8': '8', None: None}[c]


def _sym(name, latex):
    return SR.var(name, latex_name=latex)


# ---------------------------------------------------------------------------- parameters
EL, SW, CW = _sym('EL', 'e'), _sym('SW', r's_W'), _sym('CW', r'c_W')
MW, MZ, MH = _sym('MW', 'M_W'), _sym('MZ', 'M_Z'), _sym('MH', 'M_H')
GS = _sym('GS', 'g_s')
ME, MM, ML = _sym('ME', 'm_e'), _sym('MM', r'm_\mu'), _sym('ML', r'm_\tau')
MU, MC, MT = _sym('MU', 'm_u'), _sym('MC', 'm_c'), _sym('MT', 'm_t')
MD, MS, MB = _sym('MD', 'm_d'), _sym('MS', 'm_s'), _sym('MB', 'm_b')

LEPTONS = [('ne', 'e', ME), ('nm', 'mu', MM), ('nt', 'ta', ML)]
QUARKS = [('u', 'd', MU, MD), ('c', 's', MC, MS), ('t', 'b', MT, MB)]
LATEX = {'ne': r'\nu_e', 'nm': r'\nu_\mu', 'nt': r'\nu_\tau', 'e': 'e', 'mu': r'\mu', 'ta': r'\tau',
         'u': 'u', 'c': 'c', 't': 't', 'd': 'd', 's': 's', 'b': 'b'}

ALIASES = {'e-': 'e', 'e+': 'e~', 'mu-': 'mu', 'mu+': 'mu~', 'tau-': 'ta', 'tau+': 'ta~', 'tau': 'ta',
           'photon': 'A', 'gamma': 'A', 'a': 'A', 'h': 'H', 'gluon': 'g', 'nu_e': 'ne', 'nu_mu': 'nm',
           'nu_tau': 'nt', 'nu_e~': 'ne~', 'nu_mu~': 'nm~', 'nu_tau~': 'nt~',
           'ubar': 'u~', 'dbar': 'd~', 'sbar': 's~', 'cbar': 'c~', 'bbar': 'b~', 'tbar': 't~'}


def _gR(Q):
    return -SW / CW * Q


def _gL(I3, Q):
    return (I3 - SW ** 2 * Q) / (SW * CW)


def QED(muons=True, taus=True):
    """QED with e (and mu, tau): photon A, vertex -i e Q gamma^mu (Q = -1)."""
    ps = [Particle('A', 'A', 'V', 0, 0, latex=r'\gamma', xi='xiA')]
    vs = []
    for (nu, l, m) in LEPTONS[:1 + muons + taus]:
        ps.append(Particle(l, l + '~', 'F', m, -1, latex=LATEX[l]))
        vs.append(Vertex((l + '~', l, 'A'), 'FFV', I * EL, (1, 1)))
    M = Model('QED', ps, vs)
    M.aliases = dict(ALIASES)
    return M


def SM(qcd=True, xi=1):
    """The Standard Model (and QCD) of FeynArts' SM.mod / SMQCD.mod at tree level, R_xi gauge
    parameter xi for A, Z and W (xi = 1: Feynman-'t Hooft, the FeynArts default)."""
    ps = [Particle('A', 'A', 'V', 0, 0, latex=r'\gamma', xi='xiA'),
          Particle('Z', 'Z', 'V', MZ, 0, latex='Z', xi='xiZ'),
          Particle('W-', 'W+', 'V', MW, -1, latex='W^-', xi='xiW'),
          Particle('H', 'H', 'S', MH, 0, latex='H'),
          Particle('G0', 'G0', 'S', MZ, 0, latex='G^0', xi='xiZ'),
          Particle('G-', 'G+', 'S', MW, -1, latex='G^-', xi='xiW'),
          Particle('uA', 'uA~', 'U', 0, 0, latex='u_\\gamma', xi='xiA'),
          Particle('uZ', 'uZ~', 'U', MZ, 0, latex='u_Z', xi='xiZ'),
          Particle('u-', 'u-~', 'U', MW, -1, latex='u_-', xi='xiW'),
          Particle('u+', 'u+~', 'U', MW, 1, latex='u_+', xi='xiW')]
    vs = []
    # ---------------- fermions
    for (nu, l, m) in LEPTONS:
        ps += [Particle(nu, nu + '~', 'F', 0, 0, latex=LATEX[nu]), Particle(l, l + '~', 'F', m, -1, latex=LATEX[l])]
        Q, I3 = -1, -SR(1) / 2
        vs += [Vertex((l + '~', l, 'A'), 'FFV', I * EL, (-Q, -Q)),
               Vertex((nu + '~', nu, 'Z'), 'FFV', I * EL, (_gL(SR(1) / 2, 0), 0)),
               Vertex((l + '~', l, 'Z'), 'FFV', I * EL, (_gL(I3, Q), _gR(Q))),
               Vertex((nu + '~', l, 'W+'), 'FFV', I * EL / (sqrt(2) * SW), (1, 0)),
               Vertex((l + '~', nu, 'W-'), 'FFV', I * EL / (sqrt(2) * SW), (1, 0)),
               Vertex((l + '~', l, 'H'), 'FFS', -I * EL / (2 * SW * MW), (m, m)),
               Vertex((l + '~', l, 'G0'), 'FFS', -EL / (2 * SW * MW), (m, -m)),
               Vertex((nu + '~', l, 'G+'), 'FFS', -I * EL * m / (sqrt(2) * SW * MW), (0, 1)),
               Vertex((l + '~', nu, 'G-'), 'FFS', -I * EL * m / (sqrt(2) * SW * MW), (1, 0))]
    for (up, dn, mu_, md_) in QUARKS:
        ps += [Particle(up, up + '~', 'F', mu_, SR(2) / 3, '3', latex=LATEX[up]),
               Particle(dn, dn + '~', 'F', md_, -SR(1) / 3, '3', latex=LATEX[dn])]
        for (q, Q, I3, m, sgn) in ((up, SR(2) / 3, SR(1) / 2, mu_, 1), (dn, -SR(1) / 3, -SR(1) / 2, md_, -1)):
            vs += [Vertex((q + '~', q, 'A'), 'FFV', I * EL, (-Q, -Q), 'delta'),
                   Vertex((q + '~', q, 'Z'), 'FFV', I * EL, (_gL(I3, Q), _gR(Q)), 'delta'),
                   Vertex((q + '~', q, 'H'), 'FFS', -I * EL / (2 * SW * MW), (m, m), 'delta'),
                   Vertex((q + '~', q, 'G0'), 'FFS', sgn * EL / (2 * SW * MW), (m, -m), 'delta')]
            if qcd:
                vs.append(Vertex((q + '~', q, 'g'), 'FFV', -I * GS, (1, 1), 'T'))
        vs += [Vertex((up + '~', dn, 'W+'), 'FFV', I * EL / (sqrt(2) * SW), (1, 0), 'delta'),
               Vertex((dn + '~', up, 'W-'), 'FFV', I * EL / (sqrt(2) * SW), (1, 0), 'delta'),
               Vertex((up + '~', dn, 'G+'), 'FFS', I * EL / (sqrt(2) * SW * MW), (mu_, -md_), 'delta'),
               Vertex((dn + '~', up, 'G-'), 'FFS', -I * EL / (sqrt(2) * SW * MW), (md_, -mu_), 'delta')]
    # ---------------- gauge bosons
    c4 = (2, -1, -1)
    vs += [Vertex(('W+', 'W+', 'W-', 'W-'), 'VVVV', I * EL ** 2 / SW ** 2, c4),
           Vertex(('W+', 'W-', 'Z', 'Z'), 'VVVV', -I * EL ** 2 * CW ** 2 / SW ** 2, c4),
           Vertex(('W+', 'W-', 'A', 'Z'), 'VVVV', I * EL ** 2 * CW / SW, c4),
           Vertex(('W+', 'W-', 'A', 'A'), 'VVVV', -I * EL ** 2, c4),
           Vertex(('A', 'W+', 'W-'), 'VVV', -I * EL),
           Vertex(('Z', 'W+', 'W-'), 'VVV', I * EL * CW / SW)]
    # ---------------- scalars
    h4 = EL ** 2 * MH ** 2 / (SW ** 2 * MW ** 2)
    h3 = EL * MH ** 2 / (SW * MW)
    vs += [Vertex(('H', 'H', 'H', 'H'), 'SSSS', -3 * I * h4 / 4),
           Vertex(('H', 'H', 'G0', 'G0'), 'SSSS', -I * h4 / 4),
           Vertex(('H', 'H', 'G-', 'G+'), 'SSSS', -I * h4 / 4),
           Vertex(('G0', 'G0', 'G0', 'G0'), 'SSSS', -3 * I * h4 / 4),
           Vertex(('G0', 'G0', 'G-', 'G+'), 'SSSS', -I * h4 / 4),
           Vertex(('G-', 'G-', 'G+', 'G+'), 'SSSS', -I * h4 / 2),
           Vertex(('H', 'H', 'H'), 'SSS', -3 * I * h3 / 2),
           Vertex(('H', 'G0', 'G0'), 'SSS', -I * h3 / 2),
           Vertex(('G-', 'H', 'G+'), 'SSS', -I * h3 / 2)]
    e2 = EL ** 2
    vs += [Vertex(('H', 'H', 'W-', 'W+'), 'SSVV', I * e2 / (2 * SW ** 2)),
           Vertex(('G0', 'G0', 'W-', 'W+'), 'SSVV', I * e2 / (2 * SW ** 2)),
           Vertex(('G-', 'G+', 'W-', 'W+'), 'SSVV', I * e2 / (2 * SW ** 2)),
           Vertex(('G-', 'G+', 'Z', 'Z'), 'SSVV', I * e2 * (SW ** 2 - CW ** 2) ** 2 / (2 * CW ** 2 * SW ** 2)),
           Vertex(('G-', 'G+', 'A', 'Z'), 'SSVV', I * e2 * (SW ** 2 - CW ** 2) / (CW * SW)),
           Vertex(('G-', 'G+', 'A', 'A'), 'SSVV', 2 * I * e2),
           Vertex(('H', 'H', 'Z', 'Z'), 'SSVV', I * e2 / (2 * CW ** 2 * SW ** 2)),
           Vertex(('G0', 'G0', 'Z', 'Z'), 'SSVV', I * e2 / (2 * CW ** 2 * SW ** 2)),
           Vertex(('H', 'G+', 'W-', 'Z'), 'SSVV', -I * e2 / (2 * CW)),
           Vertex(('H', 'G-', 'W+', 'Z'), 'SSVV', -I * e2 / (2 * CW)),
           Vertex(('H', 'G-', 'W+', 'A'), 'SSVV', -I * e2 / (2 * SW)),
           Vertex(('H', 'G+', 'W-', 'A'), 'SSVV', -I * e2 / (2 * SW)),
           Vertex(('G-', 'G0', 'Z', 'W+'), 'SSVV', e2 / (2 * CW)),
           Vertex(('G+', 'G0', 'Z', 'W-'), 'SSVV', -e2 / (2 * CW)),
           Vertex(('G-', 'G0', 'A', 'W+'), 'SSVV', e2 / (2 * SW)),
           Vertex(('G+', 'G0', 'A', 'W-'), 'SSVV', -e2 / (2 * SW)),
           Vertex(('G0', 'H', 'Z'), 'SSV', EL / (2 * CW * SW)),
           Vertex(('G+', 'G-', 'A'), 'SSV', -I * EL),
           Vertex(('G+', 'G-', 'Z'), 'SSV', -I * EL * (SW ** 2 - CW ** 2) / (2 * CW * SW)),
           Vertex(('G-', 'H', 'W+'), 'SSV', -I * EL / (2 * SW)),
           Vertex(('G+', 'H', 'W-'), 'SSV', I * EL / (2 * SW)),
           Vertex(('G-', 'G0', 'W+'), 'SSV', EL / (2 * SW)),
           Vertex(('G+', 'G0', 'W-'), 'SSV', EL / (2 * SW)),
           Vertex(('H', 'W+', 'W-'), 'SVV', I * EL * MW / SW),
           Vertex(('H', 'Z', 'Z'), 'SVV', I * EL * MW / (SW * CW ** 2)),
           Vertex(('G+', 'W-', 'Z'), 'SVV', -I * EL * MW * SW / CW),
           Vertex(('G-', 'W+', 'Z'), 'SVV', -I * EL * MW * SW / CW),
           Vertex(('G+', 'W-', 'A'), 'SVV', -I * EL * MW),
           Vertex(('G-', 'W+', 'A'), 'SVV', -I * EL * MW)]
    # ---------------- ghosts (the gauge parameters xi are set to `xi`)
    x = SR(xi)
    vs += [Vertex(('u-~', 'u-', 'A'), 'UUV', -I * EL, (1, 0)),
           Vertex(('u+~', 'u+', 'A'), 'UUV', I * EL, (1, 0)),
           Vertex(('u-~', 'u-', 'Z'), 'UUV', I * EL * CW / SW, (1, 0)),
           Vertex(('u+~', 'u+', 'Z'), 'UUV', -I * EL * CW / SW, (1, 0)),
           Vertex(('u-~', 'uZ', 'W-'), 'UUV', -I * EL * CW / SW, (1, 0)),
           Vertex(('uZ~', 'u-', 'W+'), 'UUV', -I * EL * CW / SW, (1, 0)),
           Vertex(('u+~', 'uZ', 'W+'), 'UUV', I * EL * CW / SW, (1, 0)),
           Vertex(('uZ~', 'u+', 'W-'), 'UUV', I * EL * CW / SW, (1, 0)),
           Vertex(('u-~', 'uA', 'W-'), 'UUV', I * EL, (1, 0)),
           Vertex(('uA~', 'u-', 'W+'), 'UUV', I * EL, (1, 0)),
           Vertex(('u+~', 'uA', 'W+'), 'UUV', -I * EL, (1, 0)),
           Vertex(('uA~', 'u+', 'W-'), 'UUV', -I * EL, (1, 0)),
           Vertex(('H', 'uZ~', 'uZ'), 'SUU', -I * EL * MZ * x / (2 * SW * CW)),
           Vertex(('H', 'u-~', 'u-'), 'SUU', -I * EL * MW * x / (2 * SW)),
           Vertex(('H', 'u+~', 'u+'), 'SUU', -I * EL * MW * x / (2 * SW)),
           Vertex(('G0', 'u+~', 'u+'), 'SUU', EL * MW * x / (2 * SW)),
           Vertex(('G0', 'u-~', 'u-'), 'SUU', -EL * MW * x / (2 * SW)),
           Vertex(('G+', 'uZ~', 'u-'), 'SUU', I * EL * MZ * x / (2 * SW)),
           Vertex(('G-', 'uZ~', 'u+'), 'SUU', I * EL * MZ * x / (2 * SW)),
           Vertex(('G+', 'u+~', 'uZ'), 'SUU', I * EL * (SW ** 2 - CW ** 2) * MW * x / (2 * CW * SW)),
           Vertex(('G-', 'u-~', 'uZ'), 'SUU', I * EL * (SW ** 2 - CW ** 2) * MW * x / (2 * CW * SW)),
           Vertex(('G+', 'u+~', 'uA'), 'SUU', I * EL * MW * x),
           Vertex(('G-', 'u-~', 'uA'), 'SUU', I * EL * MW * x)]
    # ---------------- QCD
    if qcd:
        ps += [Particle('g', 'g', 'V', 0, 0, '8', latex='g', xi='xiG'),
               Particle('ug', 'ug~', 'U', 0, 0, '8', latex='u_g', xi='xiG')]
        vs += [Vertex(('g', 'g', 'g'), 'VVV', GS, (1,), 'f'),
               Vertex(('g', 'g', 'g', 'g'), 'VVVV', -I * GS ** 2, (1, 1, 1), 'ffff'),
               Vertex(('ug~', 'ug', 'g'), 'UUV', GS, (1, 0), 'f')]
    M = Model('SM' + ('QCD' if qcd else ''), ps, vs)
    M.aliases = dict(ALIASES)
    M.aliases.update({'W': 'W-', 'G': 'G-'})
    M.xi = SR(xi)
    return M
