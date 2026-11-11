"""Polar coefficients a_m(D',l') = c^F_m(n,l) + lam_m delta_m(n,l) of phi_m/Delta (Section 6), and a
check of the exact formula: its terms with c=1 against the exact coefficient at the largest n.
Usage: python3 polar.py [m ...] [--write]   (--write stores the polar terms in data/polar.json)."""
import os, sys, json, gzip
import mpmath as mp
from fractions import Fraction as Fr
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
import classnum as C
mp.mp.dps=50
def load(m):
    return {tuple(map(int,k.split(','))):v for k,v in json.load(gzip.open(os.path.join(HERE,'data',f'cF_{m}.json.gz'),'rt')).items()}
def polar(m,cF):
    out=[]
    for (n,l),v in cF.items():
        D=4*m*n-l*l
        if D<0:
            a=v+C.lam(m)*C.delta(m,n,l)
            if a: out.append((D,l,a))
    return sorted(out)
def main(m,pol,D,nu):
    k=mp.mpf(-21)/2; s=0
    for (Dp,l,a) in pol:
        for mu in ([l] if l%m==0 else [l,-l]):
            s+=mp.mpf(a.numerator)/a.denominator*mp.cos(mp.pi*mu*nu/m)*(mp.mpf(-Dp)/D)**((1-k)/2)*mp.besseli(1-k,mp.pi*mp.sqrt(-Dp*D)/m)
    return -2*mp.pi/mp.sqrt(2*m)*s
if __name__=="__main__":
    write='--write' in sys.argv
    ms=[int(a) for a in sys.argv[1:] if a!='--write'] or C.MSET
    path=os.path.join(HERE,'data','polar.json')
    stored=json.load(open(path)) if os.path.exists(path) else {}
    allok=True
    for m in ms:
        cF=load(m); pol=polar(m,cF)
        if write: stored[str(m)]=[[D,l,str(a)] for D,l,a in pol]
        same=[[D,l,str(a)] for D,l,a in pol]==stored.get(str(m))
        allok&=same
        # c=1 terms of the exact formula against the exact coefficient, at the largest n of the table
        n=max(n for n,l in cF)-1
        worst=max(abs(main(m,pol,4*m*n-l*l,l)/(cF[(n,l)]+C.lam(m)*C.delta(m,n,l))-1) for l in range(0,m+1))
        print("m",m,"polar terms",len(pol),"agree with data/polar.json:",same,"  c=1 terms / exact - 1 at n=%d:"%n,mp.nstr(worst,3),flush=True)
    if write: json.dump(stored,open(path,'w'))
    print("polar data:","PASSED" if allok else "FAILED")
