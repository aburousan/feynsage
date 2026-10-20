r"""
Running FORM (internal).

A Backend turns a FORM program into FORM's printed output.  The default, PersistentForm (see
persistent.py), keeps one FORM process alive and talks to it through pipes: a small trace takes about
0.1 ms instead of 7.6 ms with a new FORM process per call.  SubprocessBackend starts FORM for every
call.  Both have the same run() method, so the public functions do not change with the backend.
Choose with set_backend("persistent" | "subprocess") or the environment variable FEYNSAGE_BACKEND.
"""
import os, re, subprocess, tempfile, itertools
from .. import form as _form


class FormError(RuntimeError):
    """FORM stopped with an error.  .program is the FORM program that was run, .output FORM's reply."""

    def __init__(self, message, program="", output=""):
        """Store the error message with the generated FORM program and its output."""
        super().__init__(message)
        self.program, self.output = program, output


def explain_form_error(output):
    """A short Python-side explanation of a FORM error message."""
    text = output or ''
    hints = [
        (r'Illegal position for g5_|gamma5.*n dimensions|tracen', "gamma5 cannot be traced in d dimensions by FORM; use dim=4"),
        (r'Workspace overflow|WorkSpace', "FORM ran out of workspace: the expression is too large for the default settings"),
        (r'[Uu]ndefined|not declared|Unknown', "a name was not declared; please report this as a feynsage bug"),
        (r'Illegal character|syntax', "the generated FORM program has a syntax error; please report this as a feynsage bug"),
    ]
    for pat, msg in hints:
        if re.search(pat, text):
            return msg
    return "FORM stopped with an error"


class SubprocessBackend:
    """One FORM process per run (form, or tform -w threads)."""
    name = 'subprocess'

    def __init__(self):
        """Prepare a backend with a lazily created work directory and numbered jobs."""
        self._dir = None
        self._count = itertools.count()

    def _workdir(self):
        """Return the temporary job directory, creating it if needed."""
        if self._dir is None or not os.path.isdir(self._dir):
            self._dir = tempfile.mkdtemp(prefix='feynsage_form_')
        return self._dir

    def run(self, program, threads=1):
        """Run the FORM source string program and return its standard output as text.
        Use tform for threads > 1 when available, otherwise use FORM.
        Raise FormError if the executable is missing or the job fails."""
        cmd = [_form.FORM, '-q']
        if threads and threads > 1:
            tf = _form._tform()
            if tf:
                cmd = [tf, '-q', '-w%d' % threads]
        d = self._workdir()
        f = os.path.join(d, 'job%d_%d.frm' % (os.getpid(), next(self._count)))
        with open(f, 'w') as fh:
            fh.write(program)
        try:
            r = subprocess.run(cmd + [f], capture_output=True, text=True, cwd=d)
        except FileNotFoundError:
            raise FormError("FORM was not found.  Install it (see install.sh) or set FEYNSAGE_FORM to the form binary.")
        finally:
            try:
                os.remove(f)
            except OSError:
                pass
        if r.returncode != 0:
            out = r.stdout + r.stderr
            raise FormError(explain_form_error(out) + "\n--- FORM program ---\n" + program + "\n--- FORM output ---\n" + out,
                            program, out)
        return r.stdout


_BACKEND = [None]


def _default():
    """The persistent FORM process, unless FEYNSAGE_BACKEND=subprocess."""
    if os.environ.get('FEYNSAGE_BACKEND', 'persistent') == 'subprocess':
        return SubprocessBackend()
    from .persistent import PersistentForm
    return PersistentForm()


def backend():
    """Return the shared FORM backend, creating it on first use.
    FEYNSAGE_BACKEND=subprocess selects a new process per job, otherwise the default is persistent.
    A backend chosen with set_backend() is reused."""
    if _BACKEND[0] is None:
        _BACKEND[0] = _default()
    return _BACKEND[0]


def set_backend(kind="subprocess"):
    """"persistent" (one FORM process kept alive, the default) or "subprocess" (a new FORM process for
    every call).  A Backend instance with a run(program, threads) method works too."""
    if kind == "subprocess":
        _BACKEND[0] = SubprocessBackend()
    elif kind == "persistent":
        from .persistent import PersistentForm
        _BACKEND[0] = PersistentForm()
    elif hasattr(kind, 'run'):
        _BACKEND[0] = kind
    else:
        raise ValueError('backend must be "subprocess" or "persistent"')
    return _BACKEND[0]
