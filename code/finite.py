"""Finite part of Theorem B: for m in M and D > 1000 the fixed-J index d(n,J,m) is positive,
given the tail lemma (tail.py) for D >= D*(m).

For D > 1000 >= 4m^2 the first charge (n,J,m) has reduced form with c = m, so d*(n,J,m) = F_m(D).
The second charge (n,J+2,m), of discriminant D' = D - 4(J+1):
  * D' <= 0: d* = 0;
  * D' >= D*: F_m(D) > F_m(D') by the tail lemma;
  * 0 < D' < D*: either its reduced form has c = m (then d* = F_m(D')), or its reduced form
    (a,b,c) has a <= m (successive minima).  Hence d* <= M_m below.
The tail lemma gives F_m(D) > F_m(D1) whenever D >= D1 + 4 and D1 >= D*.  We find D1 in
[D*, D_end] with F_m(D1) > M_m and check directly the charges with 1000 < D < D1 + 4."""
import json, glob, os
HERE=os.path.dirname(os.path.abspath(__file__))
cF={}
for f in glob.glob(os.path.join(HERE,'data','cF_*.json')):
    m=int(f.split('_')[-1][:-5]); cF[m]={tuple(map(int,k.split(','))):v for k,v in json.load(open(f)).items()}
DSTAR={1:680,2:371,3:292,5:258,7:313,11:720,13:1030}
def reduce(a,b,c):
    while True:
        if c>a: a,c=c,a
        if abs(b)>c:
            k=(b+c)//(2*c); a,b=a-b*k+c*k*k,b-2*c*k; continue
        return a,abs(b),c
def dstar(n,l,m):
    D=4*m*n-l*l
    if D<=0: return 0
    a,b,c=reduce(n,l,m)
    return (-1)**(b+1)*cF[c][(a,b)]
def Fm(m):
    F={}
    for (n,l),v in cF[m].items():
        D=4*m*n-l*l
        if 0<=l<=m and D>0: F[D]=(-1)**(l+1)*v
    return F
allok=True
for m,Ds in DSTAR.items():
    F=Fm(m)
    assert 4*m*m<=1000
    Mred=max((-1)**(b+1)*cF[c][(a,b)] for c in range(1,m+1) for a in range(c,m+1) for b in range(0,c+1) if 4*a*c-b*b>0)
    Mlow=max([F[D] for D in F if D<Ds] or [0])
    M=max(Mred,Mlow)
    D1=min(D for D in F if D>=Ds and F[D]>M)
    # direct check of the charges with 1000 < D < D1+4 (only J with D'>0 matter; J < D/4)
    direct=0; ok=True
    for D in range(1001,D1+4):
        for J in range(0,D//4+1):
            if (D+J*J)%(4*m): continue
            n=(D+J*J)//(4*m)
            d=dstar(n,J,m)-dstar(n,J+2,m); direct+=1
            if d<=0: ok=False
    # Corollary C below D*: every reduced form [a,b,m] with 0<D<D* has (-1)^(b+1) c^F_m(a,b) > 0
    red=[(-1)**(b+1)*cF[m][(a,b)] for a in range(m,Ds//(4*m)+m+2) for b in range(0,m+1) if 0<4*a*m-b*b<Ds]
    redok=all(v>0 for v in red)
    good=ok and F[D1]>M and D1>=Ds and redok
    allok&=good
    print(f"m={m:2d} D*={Ds:4d} M_m={float(M):.4e} D1={D1:4d} F(D1)/M_m={F[D1]/M:.3f} direct charges={direct} reduced c=m forms below D*: {len(red)} positive={redok} ok={good}")
print("finite lemma:", "PASSED" if allok else "FAILED")
