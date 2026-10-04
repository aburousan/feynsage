r"""
Using every core.

Python runs one thread of Python code at a time, so work that is Python and Sage arithmetic is
spread over processes instead: the worker processes are made by fork(), so they start with
everything the main process has (the families, the reduced systems, the closures), and only the
inputs and results travel between them.  FORM, which is a separate program, can run several
threads itself (tform).

    from feynsage.parallel import pmap, cores, set_nproc
    cores()                       # the cores feynsage will use (FEYNSAGE_NPROC, else all)
    set_nproc(4)                  # use 4 from now on (1 switches parallel work off)
    pmap(f, items)                # [f(x) for x in items], in parallel

Where feynsage uses this by itself (nproc=None means "decide"):
    ibp_reduce / reduce_ff        the finite-field sample points
    quick_plot, scan              the points of a plot or a list
    Conditional.n_many            values of a symbolic C0/D0 closed form at many points
    c0_values, d0_values          many scalar triangles or boxes
    loop_many                     many one-loop integrals
    form.run_form(threads=...)    FORM with tform -w threads

A job is split only when that pays: for a few items, or items that take a millisecond each, the
cost of starting the workers (about 20 ms) is larger than the gain, and the work stays serial.
"""
import os
import time

_NPROC = None


def cores():
    """The number of worker processes feynsage uses: set_nproc(), else FEYNSAGE_NPROC, else all cores."""
    if _NPROC is not None:
        return _NPROC
    env = os.environ.get('FEYNSAGE_NPROC')
    if env:
        try:
            return max(1, int(env))
        except ValueError:
            pass
    return max(1, os.cpu_count() or 1)


def set_nproc(n):
    """Use n worker processes from now on (n = 1: everything serial; None: back to the default)."""
    global _NPROC
    _NPROC = None if n is None else max(1, int(n))


def resolve(nproc):
    """nproc as given by a caller: None or 'auto' -> cores(); a number -> that number."""
    if nproc is None or nproc == 'auto':
        return cores()
    return max(1, int(nproc))


_TASK = None


def _run(i):
    f, items = _TASK
    return f(items[i])


def pmap(f, items, nproc=None, min_items=2, probe=True, min_seconds=0.002):
    r"""
    [f(x) for x in items], spread over worker processes (fork).  f may be any function, also a
    lambda or a closure: the workers are copies of this process.  The results must be picklable
    (numbers, Sage expressions and ring elements are).

    nproc        None = cores(); 1 = serial
    probe        run f on the first item here first and go parallel only if it took longer than
                 min_seconds (then the items are many or slow enough for the workers to pay)
    Falls back to serial if the workers cannot start (no fork on the platform).
    """
    global _TASK
    items = list(items)
    n = min(resolve(nproc), len(items))
    try:
        import multiprocessing as mproc
        if mproc.current_process().daemon:          # inside a worker: workers cannot start workers
            n = 1
    except ImportError:
        n = 1
    if n <= 1 or len(items) < min_items:
        return [f(x) for x in items]
    out = []
    start = 0
    if probe:
        t0 = time.time()
        out.append(f(items[0]))
        start = 1
        if time.time() - t0 < min_seconds:
            return out + [f(x) for x in items[1:]]
    try:
        import multiprocessing as mproc
        ctx = mproc.get_context('fork')
    except (ValueError, ImportError):
        return out + [f(x) for x in items[start:]]
    _TASK = (f, items)
    try:
        with ctx.Pool(n) as pool:
            out += pool.map(_run, range(start, len(items)), chunksize=max(1, (len(items) - start) // (4 * n)))
    finally:
        _TASK = None
    return out
