# Exact decomposition Delta*psi^F_m = sum x_i (weak Jacobi monomials) - lam * G_m (Theorem D), verified on all
# coefficients n<=K, with c^F_m computed from the definition (jacobi.compute), independently of the stored tables.
# Also: the stored tables for m<=3 (computed through Theorem D) agree with the definition up to n=K2.
import os
HERE=os.path.dirname(os.path.abspath(__file__))
import json,sys,gzip
import classnum
from fractions import Fraction as Fr
from flint import fmpz_poly, fmpq_mat, fmpq
import jacobi as fast
from jacobi import B, mono_series
exec(open(os.path.join(HERE,'hurwitz.py')).read())
def sigma(n,k): return sum(d**k for d in range(1,n+1) if n%d==0)
def series(m,K):
    L=2*int((4*m*K+m*m)**0.5)+4*m+10
    E4=mono_series(K,L,{(0,0):1,**{(n,0):240*sigma(n,3) for n in range(1,K+1)}})
    E6=mono_series(K,L,{(0,0):1,**{(n,0):-504*sigma(n,5) for n in range(1,K+1)}})
    A=mono_series(K,L,{(0,1):1,(0,0):-2,(0,-1):1})
    from math import comb
    for n in range(1,K+1):
        for b in (1,-1):
            g={}
            for k in range(0,3):
                if n*k<=K: g[(n*k,b*k)]=comb(2,k)*(-1)**k
            A=A.mul(mono_series(K,L,g))
        g={};k=0
        while n*k<=K: g[(n*k,0)]=comb(k+3,3); k+=1
        A=A.mul(mono_series(K,L,g))
    c=fast.load_c()
    Bt={}
    for n in range(0,K+1):
        for l in range(-2*n-2,2*n+3):
            D=4*n-l*l
            if D>=-1 and c.get(D,0): Bt[(n,l)]=c[D]
    B2=mono_series(K,L,Bt)   # = 2 phi_{0,1}
    return E4,E6,A,B2,L
def decompose(m,K,cF,corr='G'):
    E4,E6,A,B2,L=series(m,K)
    mons=[];vecs=[]
    for i in range(0,m+1):
        j=m-i
        for a in range(0,m+2):
            for b in range(0,m+2):
                if 4*a+6*b-2*i==2:
                    f=mono_series(K,L,{(0,0):1})
                    for _ in range(a): f=f.mul(E4)
                    for _ in range(b): f=f.mul(E6)
                    for _ in range(i): f=f.mul(A)
                    for _ in range(j): f=f.mul(B2)
                    mons.append((a,b,i,j)); vecs.append(f)
    keys=[(n,l) for n in range(0,K+1) for l in range(0,m+1)]
    P=[1]+[0]*(K+3)
    for n in range(1,K+4):
        for _ in range(24):
            for i in range(n,K+4): P[i]+=P[i-n]
    Delta=[0]*(K+3)  # q*prod(1-q^n)^24 coefficients via inverse of P
    # coefficient of Delta*psiF at (n,l): sum_j tau(j) cF(n-j,l), tau(j)=coef q^j of Delta
    inv=[0]*(K+3); inv[0]=1
    for i in range(1,K+3): inv[i]=-sum(P[k]*inv[i-k] for k in range(1,i+1))
    tau=[0]+inv[:K+2]
    DF={}
    for (n,l) in keys:
        DF[(n,l)]=sum(tau[j]*cF.get((n-j,l),0) for j in range(1,n+2))
    def HVm(n,l):
        if corr=='plain':
            D=4*m*n-l*l
            if D<0: return Fr(0)
            s=Fr(0)
            for d in range(1,m+1):
                if m%d==0 and n%d==0 and l%d==0: s+=d*Hcl(D//(d*d))
            return s
        return classnum.g(m,n,l)
    cols=[[fmpq(v.coef(n,l)) for (n,l) in keys] for v in vecs]+[[fmpq(HVm(n,l).numerator,HVm(n,l).denominator) for (n,l) in keys]]
    M=fmpq_mat(len(keys),len(cols),[cols[j][i] for i in range(len(keys)) for j in range(len(cols))])
    rhs=fmpq_mat(len(keys),1,[fmpq(DF[k]) for k in keys])
    rank=M.rank()
    # solve least squares exactly via normal equations then verify
    Mt=M.transpose(); x=(Mt*M).solve(Mt*rhs)
    resid=M*x-rhs
    ok=all(resid[i,0]==0 for i in range(len(keys)))
    return mons,[x[i,0] for i in range(len(cols))],rank,len(cols),ok,len(keys)
if __name__=="__main__":
    allok=True
    for m in (list(map(int,sys.argv[1:])) or [1,2,3,4,5,7,8,9,11,13]):
        K=14
        cF=fast.compute(m,K+2)
        mons,x,rank,ncols,ok,nkeys=decompose(m,K,cF)
        lam=-x[-1]; expect=classnum.lam(m)
        good=ok and rank==ncols and lam==expect
        allok&=good
        print(f"m={m}: weak Jacobi monomials {len(mons)}, rank {rank}/{ncols}, coefficients checked {nkeys}, exact={ok}, lambda={lam} (expected {expect})")
        if m in (1,2,3,5,7,11,13):
            print("   phi_m =", " + ".join(f"({xi})*E4^{a}E6^{b}A^{i}(2phi01)^{j}" for (a,b,i,j),xi in zip(mons,x[:-1])))
    # without the correction -p H|U_p V_{p^{k-2}}, or with two prime factors, there is no such decomposition
    for m in (4,6):
        cF=fast.compute(m,14)
        mons,x,rank,ncols,ok,nkeys=decompose(m,12,cF,corr='plain')
        allok&=(not ok) and rank==ncols
        print(f"m={m}: with H|V_m alone the linear system is inconsistent (exact={ok})")
    # the stored tables for m<=3 come from Theorem D; compare with the definition
    for m,K2 in ((1,120),(2,80),(3,60)):
        st={tuple(map(int,k.split(','))):v for k,v in json.load(gzip.open(os.path.join(HERE,'data',f'cF_{m}.json.gz'),'rt')).items()}
        d=fast.compute(m,K2)
        same=all(st[k]==v for k,v in d.items() if k[0]<K2-1)
        allok&=same
        print(f"m={m}: stored table agrees with the definition for n<{K2-1}: {same}")
    # the series y S^2 used in jacobi.py equals -eta^6/theta_1^2 = prod(1-q^n)^4 / ((y^(1/2)-y^(-1/2))^2 prod(1-q^n y)^2 (1-q^n/y)^2)
    K3,W=40,40; L3=W+K3+4
    def mulp(a,b):
        out={}
        for (i,k),u in a.items():
            for (j,l),v in b.items():
                if i+j<=K3 and -L3<=k+l<=L3: out[(i+j,k+l)]=out.get((i+j,k+l),0)+u*v
        return out
    prod={(0,k+1):k+1 for k in range(0,L3)}          # y/(1-y)^2
    for n in range(1,K3+1):
        f={(n*k,0):(-1)**k*[1,4,6,4,1][k] for k in range(5) if n*k<=K3}
        g={(n*k,k):k+1 for k in range(0,K3//n+1)}       # 1/(1-q^n y)^2
        h={(n*k,-k):k+1 for k in range(0,K3//n+1)}      # 1/(1-q^n/y)^2
        prod=mulp(mulp(mulp(prod,f),g),h)
    rowsS=[[0]*(2*L3+1) for _ in range(K3+1)]
    n=0
    while n*(n+1)//2<=K3:
        e0=n*(n+1)//2; sg=-1 if n%2 else 1
        for k in range(0,L3+1):
            if e0+n*k<=K3: rowsS[e0+n*k][k+L3]+=sg
        n+=1
    n=-1
    while n*(n+1)//2<=K3:
        e0=n*(n+1)//2; sg=-1 if n%2 else 1; k=1
        while e0-n*k<=K3 and k<=L3:
            rowsS[e0-n*k][-k+L3]-=sg; k+=1
        n-=1
    Sser=B(K3,L3,[fmpz_poly(r) for r in rowsS])
    yS2=Sser.mul(Sser).mul(mono_series(K3,L3,{(0,1):1}))
    same=all(yS2.coef(i,l)==prod.get((i,l),0) for i in range(K3+1) for l in range(-W,W+1))
    allok&=same
    print(f"y S^2 = -eta^6/theta_1^2 for q-degree <= {K3}, |y-degree| <= {W}: {same}")
    print("identity:","PASSED" if allok else "FAILED")
