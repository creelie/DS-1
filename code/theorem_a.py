"""Theorem A: d*(n,l,m) > 0 for every charge with 0 < 4mn-l^2 <= Dmax, and d(n,J,m) > 0 for every
five-dimensional charge with 0 < 4mn-J^2 <= Dmax.  Usage: python3 theorem_a.py 7200 [processes]
(also writes data/margins_<Dmax>.json, the least relative margin d/d* at each D, for the figures)"""
import os, sys, json, gzip, glob
from fractions import Fraction as Fr
from multiprocessing import Pool
HERE=os.path.dirname(os.path.abspath(__file__))
def load_tables():
    cF={}
    for f in glob.glob(os.path.join(HERE,'data','cF_*.json.gz')):
        m=int(os.path.basename(f)[3:-8])
        cF[m]={tuple(map(int,k.split(','))):v for k,v in json.load(gzip.open(f,'rt')).items()}
    return cF
cF=load_tables()
def reduce(a,b,c):
    """Gauss reduction of the positive definite form [a,b,c] = a x^2 + b xy + c y^2.
    Returns (a',|b'|,c') with |b'| <= c' <= a' in the same GL2(Z) class."""
    while True:
        if c>a: a,c=c,a
        if abs(b)>c:
            k=(b+c)//(2*c)            # b-2ck lies in (-c,c]
            a,b=a-b*k+c*k*k,b-2*c*k
            continue
        return a,abs(b),c
def dstar(n,l,m):
    D=4*m*n-l*l
    if D<=0: return 0
    a,b,c=reduce(n,l,m)
    assert 4*a*c-b*b==D and b<=c<=a
    return (-1)**(b+1)*cF[c][(a,b)]
def divisors(P,spf):
    ds=[1]
    while P>1:
        p=spf[P]; e=0
        while P%p==0: P//=p; e+=1
        ds=[d*p**i for d in ds for i in range(e+1)]
    return ds
def five(args):
    Dmax,J0,step=args
    Pmax=(Dmax+(Dmax//4+1)**2)//4+1
    spf=list(range(Pmax+1))
    for i in range(2,int(Pmax**0.5)+1):
        if spf[i]==i:
            for j in range(i*i,Pmax+1,i):
                if spf[j]==j: spf[j]=i
    n5=0; bad=[]; worst=(9,(0,)); mins={}
    for J in range(J0,Dmax//4+2,step):
        for D in range(1,Dmax+1):
            if (D+J*J)%4: continue
            P=(D+J*J)//4
            for m in divisors(P,spf):
                n=P//m
                d1=dstar(n,J,m); d2=dstar(n,J+2,m)
                n5+=1
                if d1-d2<=0: bad.append((m,n,J,D,d1,d2))
                elif d2>0:
                    r=Fr(d1-d2,d1)
                    if (r,(m,n,J,D))<worst: worst=(r,(m,n,J,D))
                    if r<mins.get(D,2): mins[D]=r
    return n5,bad,worst,mins
if __name__=="__main__":
    Dmax=int(sys.argv[1]); procs=int(sys.argv[2]) if len(sys.argv)>2 else 1
    bad4=[];n4=0
    for c in range(1,int((Dmax/3)**0.5)+2):
        for b in range(0,c+1):
            for a in range(c,(Dmax+b*b)//(4*c)+1):
                D=4*a*c-b*b
                if D<=0 or D>Dmax: continue
                n4+=1
                v=dstar(a,b,c)
                if v<=0: bad4.append((a,b,c,D,v))
    print("reduced 4D forms",n4,"nonpositive",bad4[:10],flush=True)
    with Pool(procs) as pool:
        res=pool.map(five,[(Dmax,i,procs) for i in range(procs)])
    n5=sum(r[0] for r in res); bad5=sum((r[1] for r in res),[]); worst=min(r[2] for r in res)
    mins={}
    for r in res:
        for D,v in r[3].items():
            if v<mins.get(D,2): mins[D]=v
    json.dump({str(D):float(v) for D,v in sorted(mins.items())},open(os.path.join(HERE,'data',f'margins_{Dmax}.json'),'w'))
    print("5D charges (m,n,J) with J <= Dmax/4+1:",n5,"nonpositive",bad5[:10],"min relative margin",float(worst[0]),worst[1])
