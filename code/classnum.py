"""The class number term of Theorem D.  For m = 1 or a prime power p^k,
   phi_m = Delta psi^F_m + lam_m G_m  is a weak Jacobi form, where
   G_1 = H,  G_p = H|V_p,  G_{p^k} = H|V_{p^k} - p H|U_p V_{p^{k-2}}  (k >= 2),
   lam_1 = 648 = 2 p24(2) and lam_m = p24(m+1) otherwise."""
import os
from fractions import Fraction as Fr
HERE=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE,'hurwitz.py')).read())
MSET=[1,2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31]
# the discriminants D*_m from which the tail certificate holds (tail.py)
DSTAR={1:680,2:371,3:292,4:281,5:258,7:313,8:388,9:481,11:720,13:1030,16:1624,17:1857,19:2375,
       23:3630,25:4370,27:5185,29:6079,31:7052}
class _P24(list):
    """p24(n), the coefficients of 1/eta^24 = sum p24(n) q^(n-1); extended on demand through
    n p24(n) = 24 sum_{k=1}^n sigma(k) p24(n-k)."""
    def __getitem__(self,n):
        while len(self)<=n:
            N=len(self); _sig.extend(sum(d for d in range(1,k+1) if k%d==0) for k in range(len(_sig),N+1))
            self.append(24*sum(_sig[k]*list.__getitem__(self,N-k) for k in range(1,N+1))//N)
        return list.__getitem__(self,n)
_sig=[0]
P24=_P24([1])
def prime_power(m):
    """(p,k) with m = p^k, k >= 1, or None."""
    for p in range(2,m+1):
        if m%p==0:
            k=0; t=m
            while t%p==0: t//=p; k+=1
            return (p,k) if t==1 else None
    return None
def lam(m): return 648 if m==1 else P24[m+1]
def corr(m):
    """G_m as a list of (coefficient, e, k): G_m = sum coefficient * H|U_e V_k."""
    if m==1: return [(1,1,1)]
    p,k=prime_power(m)
    return [(1,1,m)] if k==1 else [(1,1,m),(-p,p,m//(p*p))]
def sigma(n): return sum(d for d in range(1,n+1) if n%d==0)
def sigstar(m):
    """sum over the terms of G_m of |coefficient| * sigma(k)."""
    return sum(abs(co)*sigma(k) for co,e,k in corr(m))
_H={}
def Hc(N):
    if N<0: return Fr(0)
    if N not in _H: _H[N]=Hcl(N)
    return _H[N]
def HUV(e,k,n,l):
    """coefficient of q^n zeta^l in H|U_e V_k (index e^2 k):
       sum_{d | (k,n,l/e)} d H((4kn-(l/e)^2)/d^2), zero unless e | l."""
    if n<0 or l%e: return Fr(0)
    l1=l//e; s=Fr(0)
    for d in range(1,k+1):
        if k%d==0 and n%d==0 and l1%d==0:
            s+=d*Hc((4*k*n-l1*l1)//(d*d))
    return s
def g(m,n,l): return sum(co*HUV(e,k,n,l) for co,e,k in corr(m))
def delta(m,n,l):
    """coefficient of q^n zeta^l in G_m/Delta."""
    return sum(P24[j]*g(m,n+1-j,l) for j in range(0,n+2))
