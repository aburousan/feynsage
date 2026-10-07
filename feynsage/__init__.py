"""feynsage: Feynman-integral tools for SageMath.

Start here:
    family(...), graph(...), diagram(name)      build integral families and Feynman graphs from strings
    ibp_reduce(fam, targets)                    IBP reduction (finite fields or exact), symmetries found automatically
    loop(numerator, props..., kin=...)          one-loop tensor integrals -> A0, B0, C0, D0 (Package-X conventions)
    A0, B0, PVB, PVC, C0, D0                    one-loop scalar and Passarino-Veltman functions
    dirac_trace(slash(p)*gamma(mu)*...)         Dirac traces, spin and colour sums in textbook notation (FORM inside)
    form.compute(...)                           the older string interface to FORM
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
from .tex import eq, notation
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
from . import qft
from .qft import (momenta, lorentz_indices, dot, comp, metric, epsilon, Dot, Comp, Metric, Delta, Epsilon,
                  set_convention, convention, spacetime_dimension,
                  gamma, slash, gamma5, PL, PR, sigma, one, u, v, ubar, vbar, DiracExpr, DiracError,
                  dirac_trace, contract, conjugate, dirac_bar, spin_sum, dummy_indices, simplify_dirac,
                  polarization, conjugate_vector, polarization_sum, to_loop, to_euclidean, to_minkowski, FormError, set_backend,
                  color_indices, quark_colors, T_color, f_color, color_delta, color_trace, color_chain, color_factor)
from .expansions import laurent
from .de import derivative, diff_reduce, differential_equation
from .kira_io import read_kira
from .counting import critical_points, master_count
from .calculus import delta_integrate, principal_value, integrate_termwise
from .topologies import topologies, Topology
from .diagrams import insert_fields, conjugate_amplitude, Diagram, DiagramList
from .models import SM, QED, Model
from .expansions import residue
from .parametric import feynman_parametrize
from .sector import sector_decompose, integrate_sectors
from . import gpl
