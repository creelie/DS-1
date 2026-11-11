"""Run every check of the paper from the stored tables in data/ (a few minutes).
Requires Python 3 with python-flint (FLINT and Arb) and mpmath.
To recompute the tables from the definitions first, run  python3 phi01.py 560  and
python3 tables.py 1000  (the second takes about ten minutes)."""
import subprocess, sys, os
HERE=os.path.dirname(os.path.abspath(__file__))
steps=[("identity.py",[],"identity: PASSED"),
       ("polar.py",[],"polar data: PASSED"),
       ("tail.py",[],"tail certificate: PASSED"),
       ("finite.py",[],"finite lemma: PASSED"),
       ("theorem_a.py",["1000"],"nonpositive []"),
       ("raw_coefficients.py",[],"bad 16")]
fail=0
for f,args,token in steps:
    out=subprocess.run([sys.executable,os.path.join(HERE,f)]+args,capture_output=True,text=True,cwd=HERE)
    ok=out.returncode==0 and token in out.stdout and (f!="theorem_a.py" or out.stdout.count(token)==2)
    fail+=not ok
    print(("PASS " if ok else "FAIL ")+f); sys.stdout.flush()
    if not ok: print(out.stdout[-2000:],out.stderr[-2000:])
print("all checks passed" if fail==0 else f"{fail} check(s) failed")
