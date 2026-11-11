"""Run every check of the paper from the stored tables in data/ (about a minute on four cores).
Requires Python 3 with python-flint (FLINT and Arb) and mpmath.
To recompute the tables first, run  python3 phi01.py 3700  and  python3 tables.py 7200
(about half an hour on four cores, with four copies of tables.py running in parallel)."""
import subprocess, sys, os
HERE=os.path.dirname(os.path.abspath(__file__))
steps=[("identity.py",[],"identity: PASSED"),
       ("polar.py",[],"polar data: PASSED"),
       ("tail.py",[],"tail certificate: PASSED"),
       ("finite.py",["7200"],"finite lemma: PASSED"),
       ("theorem_a.py",["7200",str(os.cpu_count() or 1)],"nonpositive []"),
       ("raw_coefficients.py",[],"bad 16")]
fail=0
for f,args,token in steps:
    out=subprocess.run([sys.executable,os.path.join(HERE,f)]+args,capture_output=True,text=True,cwd=HERE)
    ok=out.returncode==0 and token in out.stdout and (f!="theorem_a.py" or out.stdout.count(token)==2)
    fail+=not ok
    print(("PASS " if ok else "FAIL ")+f); sys.stdout.flush()
    if not ok: print(out.stdout[-2000:],out.stderr[-2000:])
print("all checks passed" if fail==0 else f"{fail} check(s) failed")
