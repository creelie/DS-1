# c^F_m via Theorem D: Delta psi^F_m = phi_m - lam_m H|V_m, phi_m a polynomial in E4,E6,A,B (B=2phi01)
import sys,json,os,time
from fractions import Fraction as Fr
from flint import fmpz_poly, fmpq
from jacobi import B
HERE=os.path.dirname(os.path.abspath(__file__))
c2={int(k):v for k,v in json.load(open(os.path.join(HERE,'data','c2phi01.json'))).items()}
def sigma(n,k): return sum(d**k for d in range(1,n+1) if n%d==0)
def scalar(K,L,co):
    rows=[fmpz_poly([0]*L+[co[i]]) for i in range(K+1)]
    return B(K,L,rows)
def from_dict(K,L,t):
    rows=[[0]*(2*L+1) for _ in range(K+1)]
    for (i,l),v in t.items():
        if 0<=i<=K and -L<=l<=L: rows[i][l+L]+=v
    return B(K,L,[fmpz_poly(r) for r in rows])
def series(m,K):
    L=int((4*m*K+m*m)**0.5)+2
    E4=scalar(K,L,[1]+[240*sigma(n,3) for n in range(1,K+1)])
    E6=scalar(K,L,[1]+[-504*sigma(n,5) for n in range(1,K+1)])
    # A = theta1^2/eta^6 normalised so that the q^0 row is y-2+1/y
    t={}
    R=int((2*K)**0.5)+3
    for n in range(-R,R):
        for n2 in range(-R,R):
            e=(n*n+n+n2*n2+n2)//2
            if e<=K: t[(e,n+n2+1)]=t.get((e,n+n2+1),0)+(-1)**((n+n2)%2)
    T=from_dict(K,L,t)
    # divide by prod(1-q^n)^6 : multiply by p6 coefficients
    p6=[1]+[0]*K
    for n in range(1,K+1):
        for _ in range(6):
            for i in range(n,K+1): p6[i]+=p6[i-n]
    A=T.mul(scalar(K,L,p6))
    if A.coef(0,1)!=1: A=B(K,L,[-r for r in A.rows])
    assert [A.coef(0,l) for l in (-1,0,1)]==[1,-2,1]
    Bt={}
    for n in range(0,K+1):
        lm=int((4*n+1)**0.5)
        for l in range(-lm,lm+1):
            D=4*n-l*l
            assert D<=max(c2)
            if c2.get(D,0): Bt[(n,l)]=c2[D]
    B2=from_dict(K,L,Bt)
    return E4,E6,A,B2,L
def compute(m,K,phi,lam,corr=None):
    """phi: list of ((a,b,i,j),Fraction) ; corr: list of (coef,e,k) for H|U_e V_k (default [(1,1,m)])"""
    E4,E6,A,B2,L=series(m,K)
    den=1
    for _,x in phi: den=den*x.denominator//__import__('math').gcd(den,x.denominator)
    tot=None
    pw={}
    def power(f,name,n):
        key=(name,n)
        if key in pw: return pw[key]
        if n==0: r=from_dict(K,L,{(0,0):1})
        else: r=power(f,name,n-1).mul(f)
        pw[key]=r; return r
    for (a,b,i,j),x in phi:
        f=power(E4,'E4',a).mul(power(E6,'E6',b)).mul(power(A,'A',i)).mul(power(B2,'B',j))
        num=x.numerator*(den//x.denominator)
        f=B(K,L,[num*r for r in f.rows])
        tot=f if tot is None else tot.add(f)
    Hmax=4*m*K+m*m+1
    exec(open(os.path.join(HERE,'hurwitz.py')).read(),globals())
    Hc={}
    def H(N):
        if N<0: return Fr(0)
        if N not in Hc: Hc[N]=Hcl(N)
        return Hc[N]
    if corr is None: corr=[(1,1,m)]
    def HUV(e,k,n,l):
        s=Fr(0)
        for d in range(1,k+1):
            if k%d==0 and n%d==0 and l%d==0:
                n2,l2=n*k//(d*d),l//d
                if l2%e: continue
                l3=l2//e; s+=d*H(4*n2-l3*l3)
        return s
    P=[1]+[0]*(K+2)
    for n in range(1,K+3):
        for _ in range(24):
            for i in range(n,K+3): P[i]+=P[i-n]
    res={}
    for l in range(0,m+1):
        g=[Fr(tot.coef(i,l),den)-lam*sum(co*HUV(e,k,i,l) for co,e,k in corr) for i in range(K+1)]
        for n in range(-1,K):
            v=sum(P[j]*g[n+1-j] for j in range(0,n+2))
            assert v.denominator==1,(n,l,v)
            res[(n,l)]=int(v)
    return res
PHI={1:[((1,0,1,0),Fr(3))],2:[((1,0,1,1),Fr(11,3)),((0,1,2,0),Fr(-10,3))],
     3:[((1,0,1,2),Fr(467,192)),((2,0,3,0),Fr(203,48)),((0,1,2,1),Fr(-215,48))]}
LAM={1:648,2:3200,3:25650}
if __name__=="__main__":
    D0=int(sys.argv[1])
    for m in map(int,sys.argv[2:]):
        K=(D0+m*m)//(4*m)+2; t=time.time()
        r=compute(m,K,PHI[m],LAM[m])
        out=os.path.join(HERE,'data',f'cF_{m}.json')
        json.dump({f"{n},{l}":v for (n,l),v in r.items()},open(out+'.tmp','w')); os.rename(out+'.tmp',out)
        print(m,K,round(time.time()-t,1),flush=True)
