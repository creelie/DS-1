"""Finite part of Theorem B (Section 9).  For m in M let
   M_m = max( (-1)^(b+1) c^F_c(a,b) over reduced forms with c <= a <= m,
              (-1)^(l+1) c^F_m(n,l) over 0 <= l <= m, 0 < 4mn - l^2 < D*_m ),
and let D_1 >= D*_m be the least discriminant of index m at which some class has
(-1)^(l+1) c^F_m(n,l) > M_m.  The proof needs D_1 + 4 <= Dmax + 1, where Theorem A covers D <= Dmax."""
import os, sys, json, gzip, glob
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
from classnum import MSET, DSTAR as C_DSTAR
DMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7200
cF={}
for f in glob.glob(os.path.join(HERE,'data','cF_*.json.gz')):
    m=int(os.path.basename(f)[3:-8]); cF[m]={tuple(map(int,k.split(','))):v for k,v in json.load(gzip.open(f,'rt')).items()}
DSTAR=C_DSTAR
allok=True
for m in MSET:
    Ds=DSTAR[m]
    vals={}   # D -> list of (-1)^(l+1) c^F_m(n,l) over the classes 0<=l<=m with discriminant D
    for (n,l),v in cF[m].items():
        D=4*m*n-l*l
        if 0<=l<=m and D>0: vals.setdefault(D,set()).add((-1)**(l+1)*v)
    Mred=max((-1)**(b+1)*cF[c][(a,b)] for c in range(1,m+1) for a in range(c,m+1) for b in range(0,c+1) if 4*a*c-b*b>0)
    Mlow=max(max(vals[D]) for D in vals if D<Ds)
    M=max(Mred,Mlow)
    D1=min(D for D in vals if D>=Ds and max(vals[D])>M)
    good=D1>=Ds and D1+4<=DMAX+1 and 4*m*m<Ds and 4*m*m<=DMAX
    allok&=good
    print(f"m={m:2d} D*={Ds:4d} classes per D<={max(len(s) for s in vals.values())} M_m={float(M):.4e} D1={D1:4d} max F(D1)/M_m={max(vals[D1])/M:.3f} ok={good}")
print("finite lemma:", "PASSED" if allok else "FAILED")
