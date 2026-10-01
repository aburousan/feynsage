# Run a command, write its wall time, CPU time and peak memory (of the child processes) as JSON.
#     python3 measure.py out.json timeout_seconds cmd args...
import sys, time, json, resource, subprocess
out, limit, cmd = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
t0 = time.time()
try:
    rc = subprocess.call(cmd, timeout=limit)
except subprocess.TimeoutExpired:
    rc = "timeout"
wall = time.time() - t0
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
json.dump({"wall": wall, "cpu": ru.ru_utime + ru.ru_stime, "maxrss_mb": ru.ru_maxrss / 1024.0, "rc": rc},
          open(out, "w"))
