r"""
QFT algebra in Python and Sage notation, with FORM doing the work (no FORM code to write).

    from feynsage import *
    p, q = momenta("p q"); mu, nu = lorentz_indices("mu nu")
    dirac_trace(slash(p) * gamma(mu) * slash(q) * gamma(nu))

See docs/dirac_algebra.md.  The lower-level string interface form.compute stays available.
"""
from .tensors import (momenta, lorentz_indices, dot, comp, metric, epsilon, Dot, Comp, Metric, Delta, Epsilon,
                      set_convention, convention,
                      dimension, is_momentum, is_index)
from .tensors import dimension as spacetime_dimension, polarization, conjugate_vector, to_euclidean, to_minkowski
from .operations import polarization_sum
from .objects import (gamma, slash, gamma5, PL, PR, sigma, one, u, v, ubar, vbar, DiracExpr, DiracError)
from .operations import dirac_trace, contract, conjugate, dirac_bar, spin_sum, dummy_indices
from .simplify import simplify_dirac
from .backend import FormError, set_backend
from .compiler import CompileError
from ..form import to_loop
from .color import color_indices, quark_colors, T_color, f_color, color_delta, color_trace, color_chain, color_factor
