"""Recompute the tables data/cF_m.json.gz (exact coefficients c^F_m(n,l), 0<=l<=m, -1<=n<K, of psi^F_m)
for every m <= sqrt(D0/3)+1, which includes every index of a reduced form of discriminant at most D0.
For m <= 3 the coefficients are obtained from Theorem D (viaD.py), for m >= 4 from the definition
(jacobi.py); identity.py checks the first route against the second.
Usage: python3 phi01.py 3700; python3 tables.py 7200   (about half an hour on four cores with four copies running; needs python-flint).
Several copies may run at once: each index is claimed through a lock file."""
import sys, json, os, time, gzip
import jacobi, viaD
HERE=os.path.dirname(os.path.abspath(__file__))
D0=int(sys.argv[1])
M=int((D0/3)**0.5)+1
for m in range(1,M+1):
    out=os.path.join(HERE,'data',f'cF_{m}.json.gz')
    if os.path.exists(out): continue
    try: fd=os.open(out+'.lock',os.O_CREAT|os.O_EXCL)
    except FileExistsError: continue
    K=(D0+m*m)//(4*m)+2
    t=time.time()
    r=viaD.compute(m,K,viaD.PHI[m],viaD.LAM[m]) if m<=3 else jacobi.compute(m,K)
    with gzip.open(out+'.tmp','wt') as f: json.dump({f"{n},{l}":v for (n,l),v in r.items()},f)
    os.rename(out+'.tmp',out); os.close(fd); os.remove(out+'.lock')
    print(m,K,round(time.time()-t,1),flush=True)
