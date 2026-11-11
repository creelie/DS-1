"""Exact coefficients c^F_m(n,l), 0<=l<=m, -1<=n<=K, of the finite part psi^F_m of the
m-th Fourier-Jacobi coefficient of 1/Phi_10.
Method: Delta*psi_m = chi(Sym^{m+1}K3) * Kf with Kf = -eta^6/theta_1^2 expanded in |q|<|y|<1;
psi^F_m = psi_m - p24(m+1) A_{1,m}/Delta.  Bivariate series = list of rows (fmpz_poly in y,
offset L), row i = coefficient of q^i."""
import sys, json, os
from flint import fmpz_poly, fmpz
from math import comb
def load_c():
    return {int(k):v for k,v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','c2phi01.json'))).items()}
class B:
    def __init__(s,K,L,rows=None):
        s.K=K; s.L=L; s.rows=rows if rows is not None else [fmpz_poly(0) for _ in range(K+1)]
    def trunc_row(s,p):
        # p has offset 2L after product; bring back to offset L, keep |l|<=L
        c=p.coeffs(); L=s.L
        return fmpz_poly(c[L:3*L+1])
    def mul(s,o):
        # Kronecker substitution: row i, offset index t in [0,2L] -> x^(i*W+t), W=4L+1
        K=s.K; L=s.L; W=4*L+1
        def pack(rows):
            c=[]
            for i in range(K+1):
                r=[int(x) for x in rows[i].coeffs()]
                r=r+[0]*(W-len(r))
                c.extend(r)
            return fmpz_poly(c)
        P=pack(s.rows)*pack(o.rows)
        c=P.coeffs(); out=[]
        for i in range(K+1):
            seg=c[i*W+L:i*W+3*L+1]
            out.append(fmpz_poly(seg))
        return B(K,L,out)
    def add(s,o,f=1):
        return B(s.K,s.L,[a+f*b for a,b in zip(s.rows,o.rows)])
    def scal_div(s,d):
        out=[]
        for a in s.rows:
            cs=[int(x) for x in a.coeffs()]
            assert all(x%d==0 for x in cs)
            out.append(fmpz_poly([x//d for x in cs]))
        return B(s.K,s.L,out)
    def coef(s,i,l):
        c=s.rows[i].coeffs(); j=l+s.L
        return int(c[j]) if 0<=j<len(c) else 0
def mono_series(K,L,terms):
    rows=[[0]*(2*L+1) for _ in range(K+1)]
    for (i,l),v in terms.items():
        if 0<=i<=K and -L<=l<=L: rows[i][l+L]+=v
    return B(K,L,[fmpz_poly(r) for r in rows])
def compute(m,K,L=None):
    c=load_c()
    Nmax=m+1
    if L is None: L=int(2*((Nmax)*(K+1))**0.5)+m+K+8
    def cc(D):
        if D< -1: return 0
        assert D<=max(c),("need c up to",D)
        return c.get(D,0)
    # kL_k series: coefficient at (n,l) = sum_{r|gcd(k,n,l)} (k/r) c(4nk/r^2 - (l/r)^2)
    from math import gcd
    kLk={}
    for k in range(1,Nmax+1):
        t={}
        for n in range(0,K+1):
            lm=int((4*n*k+k*k)**0.5)+1
            for l in range(-lm,lm+1):
                g=gcd(gcd(k,n),abs(l)) if (n or l) else k
                s=0
                for r in range(1,g+1):
                    if k%r==0 and n%r==0 and l%r==0:
                        s+=(k//r)*cc(4*(n//r)*(k//r)-(l//r)**2)
                if s: t[(n,l)]=s
        kLk[k]=mono_series(K,L,t)
    Z=[mono_series(K,L,{(0,0):1})]
    for N in range(1,Nmax+1):
        acc=B(K,L)
        for k in range(1,N+1):
            acc=acc.add(kLk[k].mul(Z[N-k]))
        Z.append(acc.scal_div(N))
    chi=Z[Nmax]
    # Kf = -eta^6/theta_1^2 = y*S^2,  S = sum_n (-1)^n q^{n(n+1)/2}/(1-q^n y)  (|q|<|y|<1)
    rows=[[0]*(2*L+1) for _ in range(K+1)]
    n=0
    while n*(n+1)//2<=K:
        e0=n*(n+1)//2; sg=(-1 if n%2 else 1)
        if n==0:
            for k in range(0,L+1): rows[0][k+L]+=1
        else:
            k=0
            while e0+n*k<=K:
                if k<=L: rows[e0+n*k][k+L]+=sg
                k+=1
        n+=1
    n=-1
    while n*(n+1)//2<=K:
        e0=n*(n+1)//2; sg=(-1 if n%2 else 1); k=1
        while e0-n*k<=K:
            if -k>=-L: rows[e0-n*k][-k+L]-=sg
            k+=1
        n-=1
    Sser=B(K,L,[fmpz_poly(r) for r in rows])
    Kf=Sser.mul(Sser).mul(mono_series(K,L,{(0,1):1}))
    DP=chi.mul(Kf)   # Delta * psi_m, coefficient (i,l) <-> q^{i} y^l times q^0 ... Delta=q*prod => Delta*psi_m has q^{>=0}
    # p24
    P=[1]+[0]*(K+2)
    for n in range(1,K+3):
        for _ in range(24):
            for i in range(K+2,n-1,-1): P[i]+=P[i-n]
    # wait: that's prod (1-q^n)^{-1} applied 24 times: P[i]+=P[i-n] ascending order needed for 1/(1-q^n)
    KP=max(K,m)+3
    P=[1]+[0]*(KP-1)
    for n in range(1,KP):
        for _ in range(24):
            for i in range(n,KP): P[i]+=P[i-n]
    p24=P
    res={}
    for l in range(0,m+1):
        # g(i) = coefficient q^i y^l of Delta*psi^F_m = DP - p24(m+1)*A_{1,m}
        g=[DP.coef(i,l) for i in range(K+1)]
        # A_{1,m} coefficients at y^l, q^i (annulus expansion) : all have D<0
        for s in range(-K-2,K+3):
            e0=m*s*s+s
            if s>=0:
                k=l-(2*m*s+1)
                if k>=0 and e0+s*k<=K: g[e0+s*k]-=p24[m+1]*(k+1)
            else:
                k=2*m*s+1-2-l
                if k>=0:
                    e=e0-2*s-s*k
                    if 0<=e<=K: g[e]-=p24[m+1]*(k+1)
        # psi^F = g / Delta = q^{-1} * g * P  -> c^F(n,l) = sum_j P[j] g[n+1-j]
        for n in range(-1,K):
            res[(n,l)]=sum(P[j]*g[n+1-j] for j in range(0,n+2))
    return res
if __name__=="__main__":
    m=int(sys.argv[1]); K=int(sys.argv[2])
    r=compute(m,K)
    json.dump({f"{n},{l}":v for (n,l),v in r.items()},open(f'fcF_{m}.json','w'))
    print(m,K,len(r))
