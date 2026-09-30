"""feynsage: Feynman-integral tools for SageMath.

Start here:
    family(...), graph(...), diagram(name)      build integral families and Feynman graphs from strings
    ibp_reduce(fam, targets)                    IBP reduction (finite fields or exact), symmetries found automatically
    loop(numerator, props..., kin=...)          one-loop tensor integrals -> A0, B0, C0, D0 (Package-X conventions)
    A0, B0, PVB, PVC, C0, D0                    one-loop scalar and Passarino-Veltman functions
    form.dirac_trace(...)                       Dirac traces with FORM (gamma_5 included)
    plotting.quick_plot, g.plot()               plots and diagrams in the style of the notes
Every public function takes explain=True, or call feynsage.info(function), to print what its output means.
"""
from .momenta import mom, Kinematics
from .family import IntegralFamily
from .laporta import Reducer
from .graph import FeynmanGraph
from .ff import reduce_ff, IBPSystem
from .easy import family, graph, diagram, symmetries, ibp_reduce, kinematics
from .pv import loop, A0, B0, PVB, PVC, C0, D0, DiscB, LogM, uv_part, finite_part, c0_numeric, d0_numeric, eps, mu
from . import oneloop, form, plotting, pv, easy
from .explain import info, _install
_install()
from .pv import A0, B0, uv_part, finite_part, c0_numeric, d0_numeric
from .ff import reduce_ff
