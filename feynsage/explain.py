r"""
info(f): print what a feynsage function computes and what its output means.

Every public function also accepts explain=True, which prints the same text before
computing.  The text is the function's docstring; it always says what comes back.
"""


def info(obj):
    """Print the description of a feynsage function, class or result object."""
    if hasattr(obj, 'info') and callable(obj.info) and not isinstance(obj, type) and not callable(getattr(obj, '__wrapped__', None)):
        try:
            return obj.info()
        except TypeError:
            pass
    doc = getattr(obj, '__doc__', None) or "no description"
    name = getattr(obj, '__name__', type(obj).__name__)
    print("%s\n%s\n%s" % (name, "-" * len(name), doc.strip("\n")))


def explainable(f):
    """Give f an explain=False keyword: explain=True prints what f computes and returns, then runs it."""
    import functools
    import inspect
    try:
        if 'explain' in inspect.signature(f).parameters:
            return f
    except (TypeError, ValueError):
        pass

    @functools.wraps(f)
    def wrapper(*args, explain=False, **kw):
        if explain:
            info(f)
        return f(*args, **kw)
    return wrapper


def _install():
    """Add explain= to the public functions and methods of feynsage."""
    from . import oneloop, form, ff, pv, plotting
    from .graph import FeynmanGraph
    from .family import IntegralFamily
    from .laporta import Reducer
    for mod, names in ((oneloop, ['tadpole', 'G', 'bubble_massless', 'triangle_onshell', 'bubble_equal_mass',
                                  'bubble_one_mass', 'expand_eps']),
                       (form, ['trace', 'dirac_trace', 'to_loop', 'run_form']),
                       (ff, ['reduce_ff']),
                       (pv, ['A0', 'B0', 'C0', 'D0', 'uv_part', 'pole_parts', 'finite_part', 'c0_numeric', 'd0_numeric', 'explicit']),
                       (plotting, ['set_theme', 'curves', 'mb_plane', 'draw_graph', 'draw_panels', 'draw_sectors'])):
        for n in names:
            if hasattr(mod, n):
                setattr(mod, n, explainable(getattr(mod, n)))
    for cls, names in ((FeynmanGraph, ['U', 'F', 'F0', 'U_kirchhoff', 'spanning_trees', 'two_forests',
                                             'family', 'contract', 'plot']),
                       (IntegralFamily, ['UF', 'is_zero_sector', 'ibp']),
                       (Reducer, ['run', 'reduce', 'masters'])):
        for n in names:
            if hasattr(cls, n):
                setattr(cls, n, explainable(getattr(cls, n)))
