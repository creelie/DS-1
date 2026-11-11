"""Recompute the tables data/cF_m.json (exact coefficients of psi^F_m) for every m with
3m^2 < D0, i.e. every index that occurs in a reduced form of discriminant below D0.
Usage: python3 tables.py 1000   (about ten minutes; needs python-flint)."""
import sys, json, os, time
import jacobi
HERE=os.path.dirname(os.path.abspath(__file__))
D0=int(sys.argv[1])
M=int((D0/3)**0.5)+1
for m in range(1,M+1):
    K=(D0+m*m)//(4*m)+2
    t=time.time()
    r=jacobi.compute(m,K)
    json.dump({f"{n},{l}":v for (n,l),v in r.items()},open(os.path.join(HERE,'data',f'cF_{m}.json'),'w'))
    print(m,K,round(time.time()-t,1),flush=True)
