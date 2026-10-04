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
from .pv import loop, A0, B0, PVB, PVC, PVD, C0, D0, DiscB, LogM, uv_part, pole_parts, finite_part, c0_numeric, d0_numeric, eps, mu, explicit
from .scalar import d0_closed, d0_value
from . import ir
from .tex import eq
from .scalar import c0_closed, c0_value
from . import oneloop, form, plotting, pv, easy
from .explain import info, _install
_install()
from .pv import A0, B0, C0, D0, uv_part, pole_parts, finite_part, c0_numeric, d0_numeric, explicit
from .pv import PVA, loop_diff, loop_series, kallen, kibble, mandelstam, disc_expand, Ln, DiLog, continued_dilog
from . import dirac
from .dirac import (line_expand, loop_line, loop_matrix, projector, form_factor, transverse, longitudinal,
                    chisholm, to_chiral, to_g5, line_product)
from . import closed
from .closed import c0_expand, d0_expand, expand_c0d0, Conditional
from . import parallel
from .parallel import pmap, cores, set_nproc
from .pv import loop_many
from .scalar import c0_values, d0_values
from .plotting import scan
from .ff import reduce_ff
