r"""
A FORM process kept alive between calls (internal; choose it with set_backend("persistent")).

FORM is started as `form -pipe r,w init.frm` (the protocol of FORM's external channels, the same one
python-form uses): FORM writes its process id, we answer "pid,ourpid", init.frm sends "OK" and then
loops over #fromexternal, which reads FORM statements from our pipe up to a prompt line.  Each job is
the program feynsage would run with `form job.frm`, rewritten so that its expressions are sent back with
#toexternal instead of being printed, and dropped afterwards.  If FORM stops with an error the process
is restarted on the next call; with threads > 1 the job goes to tform in its own process.
"""
import os, re, select, subprocess, tempfile, atexit
from .. import form as _form
from .backend import FormError, SubprocessBackend, explain_form_error

_INIT = """Off statistics;
Off finalstats;
#ifndef `PIPES_'
  #message "feynsage: no pipes"
  .end
#endif
#setexternal `PIPE1_'
#toexternal "OK"
#do FSLOOP=1,1
  #fromexternal
#enddo
.end
"""
_PROMPT = 'FSREADY'
_END = 'FSJOBEND'


class PersistentForm:
    name = 'persistent'

    def __init__(self):
        self.proc = None
        self._broken = False
        self._sub = SubprocessBackend()
        self._dir = tempfile.mkdtemp(prefix='feynsage_pform_')
        self._init = os.path.join(self._dir, 'init.frm')
        with open(self._init, 'w') as fh:
            fh.write(_INIT)
        atexit.register(self.close)

    # -- process
    def _start(self):
        r_child, w_parent = os.pipe()
        r_parent, w_child = os.pipe()
        self.proc = subprocess.Popen([_form.FORM, '-q', '-pipe', '%d,%d' % (r_child, w_child), self._init],
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                     pass_fds=(r_child, w_child), cwd=self._dir)
        os.close(r_child)
        os.close(w_child)
        self.r, self.w = r_parent, w_parent
        self.owner = os.getpid()
        self.buf = ''
        line = self._read_until('\n', timeout=20)
        formpid = int(line.strip())
        os.write(self.w, ('%d,%d\n' % (formpid, os.getpid())).encode())
        self._read_until('OK', timeout=20)
        os.write(self.w, ('#prompt %s\n#-\n' % _PROMPT).encode())

    def _alive(self):
        if self.proc is not None and getattr(self, 'owner', None) != os.getpid():
            self.proc = None             # a forked child (parallel.pmap): never share the parent's pipes
        return self.proc is not None and self.proc.poll() is None

    def close(self):
        if self.proc is not None and getattr(self, 'owner', None) != os.getpid():
            self.proc = None
            return
        if self.proc is not None:
            try:
                if self.proc.poll() is None:
                    self.proc.kill()
                self.proc.wait(timeout=5)
            except Exception:
                pass
            for fd in (getattr(self, 'r', None), getattr(self, 'w', None)):
                try:
                    os.close(fd)
                except Exception:
                    pass
            self.proc = None

    def _read_until(self, mark, timeout=None):
        """Read from FORM's channel until mark; returns the text before it.  Watches FORM's stdout
        for error messages (FORM stops after an error)."""
        out_fd = self.proc.stdout.fileno()
        log = ''
        while mark not in self.buf:
            ready, _, _ = select.select([self.r, out_fd], [], [], timeout)
            if not ready:
                raise FormError("FORM did not answer in time")
            if self.r in ready:
                chunk = os.read(self.r, 1 << 16)
                if chunk:
                    self.buf += chunk.decode()
            if out_fd in ready:
                chunk = os.read(out_fd, 1 << 16)
                if chunk:
                    log += chunk.decode()
                if not chunk or '-->' in log or '==>' in log or self.proc.poll() is not None:
                    rest = ''
                    try:
                        rest = self.proc.stdout.read() or b''
                        rest = rest.decode() if isinstance(rest, bytes) else rest
                    except Exception:
                        pass
                    self.close()
                    raise FormError(explain_form_error(log + rest), "", log + rest)
        i = self.buf.index(mark)
        text, self.buf = self.buf[:i], self.buf[i + len(mark):]
        return text

    # -- jobs
    @staticmethod
    def _job(program):
        """A file program -> statements for #fromexternal that send each Local back and drop it."""
        names = re.findall(r'(?m)^\s*Local\s+(\w+)\s*=', program)
        body = re.sub(r'(?m)^\s*(Print[^;]*;|\.end)\s*$', '', program)
        if not re.search(r'(?m)^\s*Dimension\s', body):
            body = 'Dimension 4;\n' + body           # a d-dimensional job before must not leak
        lines = [body.rstrip(), '.sort']
        for n in names:
            lines.append('#toexternal "\\n%s=%%E;\\n", %s' % (n, n))
        lines += ['#toexternal "%s"' % _END, 'Drop;', '.sort', '#redefine FSLOOP "0"', _PROMPT, '']
        return '\n'.join(lines)

    def run(self, program, threads=1):
        if threads and threads > 1:
            return self._sub.run(program, threads)
        if self._broken:
            return self._sub.run(program, threads)
        job = self._job(program)
        for attempt in (0, 1):
            if not self._alive():
                try:
                    self._start()
                except (OSError, ValueError, FormError):
                    self.close()
                    self._broken = True                  # FORM without -pipe: one process per call
                    return self._sub.run(program, threads)
            try:
                os.write(self.w, job.encode())
                return self._read_until(_END)
            except (FormError, OSError) as e:
                self.close()             # a fresh FORM forgets the declarations of earlier jobs
                if attempt == 1:
                    out = getattr(e, 'output', '') or str(e)
                    raise FormError(explain_form_error(out) + "\n--- FORM program ---\n" + program
                                    + "\n--- FORM output ---\n" + out, program, out)
