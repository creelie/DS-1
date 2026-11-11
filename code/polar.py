import os
HERE=os.path.dirname(os.path.abspath(__file__))
import json, sys, mpmath as mp
from fractions import Fraction as Fr
exec(open(os.path.join(HERE,'hurwitz.py')).read())
mp.mp.dps=50
P24=[1]+[0]*400
for n in range(1,401):
    for _ in range(24):
        for i in range(n,401): P24[i]+=P24[i-n]
def LAM(m): return 648 if m==1 else P24[m+1]
def load(m):
    c=json.load(open(os.path.join(HERE,'data',f'cF_{m}.json'))); return {tuple(map(int,k.split(','))):v for k,v in c.items()}
def HVm(m,n,l):
    D=4*m*n-l*l
    if D<0 or n<0: return Fr(0)
    s=Fr(0)
    for d in range(1,m+1):
        if m%d==0 and n%d==0 and l%d==0: s+=d*Hcl(D//(d*d))
    return s
def hdelta(m,n,l): return sum(P24[j]*HVm(m,n+1-j,l) for j in range(0,n+2))
def polar(m):
    cF=load(m); lam=LAM(m)
    out=[]
    for (n,l),v in cF.items():
        D=4*m*n-l*l
        if D<0:
            a=v+lam*hdelta(m,n,l)
            if a: out.append((D,l,a))
    return sorted(out)
def main(m,pol,D,nu):
    k=mp.mpf(-21)/2; s=0
    for (Dp,l,a) in pol:
        mus=[l] if l%m==0 else [l,-l]
        for mu in mus:
            s+=mp.mpf(a.numerator)/a.denominator*mp.cos(mp.pi*mu*nu/m)*(mp.mpf(-Dp)/D)**((1-k)/2)*mp.besseli(1-k,mp.pi*mp.sqrt(-Dp*D)/m)
    return -2*mp.pi/mp.sqrt(2*m)*s
if __name__=="__main__":
    ms=list(map(int,sys.argv[1:])) or [1,2,3,5,7,11,13]
    stored=json.load(open(os.path.join(HERE,'data','polar.json')))
    allok=True
    for m in ms:
        pol=polar(m); cF=load(m); lam=LAM(m)
        same=[(D,l,str(a)) for D,l,a in pol]==[tuple(t) for t in stored[str(m)]]
        allok&=same
        print("m",m,"polar terms",len(pol),"agree with data/polar.json:",same)
        # the term c=1 of the Rademacher series against the exact coefficient (largest n in the table)
        K=max(n for n,l in cF)
        worst=max(abs(main(m,pol,4*m*(K-1)-l*l,l)/(cF[(K-1,l)]+lam*hdelta(m,K-1,l))-1) for l in range(0,m+1))
        print("   c=1 Rademacher term / exact coefficient - 1, worst over l at n=%d:"%(K-1), mp.nstr(worst,3))
    print("polar data:","PASSED" if allok else "FAILED")
