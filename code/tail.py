# Rigorous (Arb) tail certificate.  For D >= Dstar (D in a fixed residue class mod 4m):
#   |F_m(D) - T(D)| <= eps(D) T(D),  eps < 1/2, and eps(D)+eps(D+4) < 1 - T(D)/T(D+4).
import json, sys, os
from flint import arb, ctx
from fractions import Fraction as Fr
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from classnum import lam, sigstar as sig, MSET, DSTAR as C_DSTAR
ctx.prec=200
pi=arb.pi(); nu=arb(23)/2
c0=arb(1)
for n in range(1,60): c0*= (1-(-2*pi*n).exp())**(-24)
c0=c0*(1+arb('1e-100'))
Z=arb(0)
for c in range(2,400): Z+= arb(c)**(-nu)
Z+= arb(399)**(1-nu)/(nu-1)
S=arb(1)/5
def A(x): return arb(x.numerator)/x.denominator
def small(m,Dt,Dstar):
    # polar terms whose Bessel argument at Dstar is below 12 use the increasing bound B2
    return bool(pi*(arb(-Dt)*Dstar).sqrt()/m < 12)
def terms(m,pol,D,Dstar=None):
    if Dstar is None: Dstar=D
    D=arb(D)
    Dd=-m*(m+4)
    dom=[t for t in pol if t[0]==Dd]; assert len(dom)==1 and Fr(dom[0][2])==m+2 and dom[0][1]==m
    xd=pi*(arb(-Dd)*D).sqrt()/m
    kap=(1-S/2)**11*(S*xd).gamma_lower(arb(12),regularized=1)
    TL=2*pi/arb(2*m).sqrt()*(m+2)*(arb(-Dd)/D)**(nu/2)*xd.exp()/(2*pi*xd).sqrt()*kap
    E1=arb(0);E2=arb(0)
    for (Dt,l,a) in pol:
        w=1 if l%m==0 else 2
        a=A(abs(Fr(a)))
        xt=pi*(arb(-Dt)*D).sqrt()/m
        if Dt!=Dd:
            if small(m,Dt,Dstar): E1+=2*pi/arb(2*m).sqrt()*a*w*(arb(-Dt)/D)**(nu/2)*(xt/2)**nu*xt.exp()/(nu+1).gamma()
            else: E1+=2*pi/arb(2*m).sqrt()*a*w*(arb(-Dt)/D)**(nu/2)*xt.exp()/(2*pi*xt).sqrt()
        E2+=2*pi*a*w*(arb(-Dt)/D)**(nu/2)*(xt/2)**nu*(xt/2).exp()/(nu+1).gamma()*Z
    J=(D+4*m)/(4*m)
    E3=lam(m)*sig(m)*(D+4*m+1)*(25+(J-1)*c0*(J-1)**(-6)*(4*pi*(J-1).sqrt()).exp())
    alpha=pi*(arb(m+4)/m).sqrt()
    u0=alpha*D.sqrt()
    beta=u0/(nu+1+(u0*u0+(nu+1)**2).sqrt())
    gap=1-(-2*alpha*beta/(D+4).sqrt()).exp()   # lower bound for 1-T(D)/T(D+4)
    return TL,E1,E2,E3,gap
def mono(m,pol,D):
    # sqrt(D) must exceed 2a/c for every term (see paper, Lemma on monotonicity)
    Dd=-m*(m+4); ad=pi*arb(-Dd).sqrt()/m; need=[]
    for (Dt,l,a) in pol:
        at=pi*arb(-Dt).sqrt()/m
        if Dt!=Dd: need.append((13 if small(m,Dt,D) else 1)/(ad-at))   # 2a/c, a=13/2 (bound B2) or 1/2 (B1)
        need.append(13/(ad-at/2))                    # a=13/2
    need.append(15/ad)                               # a=15/2
    need.append(5/(ad-2*pi/arb(m).sqrt()))           # a=5/2
    mx=max(need, key=lambda x: float(x.mid()))
    return bool(arb(D).sqrt()>mx), mx
def ok(m,pol,D):
    TL,E1,E2,E3,gap=terms(m,pol,D)
    eps=(E1+E2+E3)/TL
    return bool(2*eps<gap) and bool(eps<arb(1)/2) and mono(m,pol,D)[0], eps, gap
if __name__=="__main__":
    import os
    HERE=os.path.dirname(os.path.abspath(__file__))
    pols=json.load(open(os.path.join(HERE,'data','polar.json')))
    DSTAR=C_DSTAR
    assert sorted(DSTAR)==MSET
    allok=True
    for m,Ds in DSTAR.items():
        pol=[(D,l,a) for D,l,a in pols[str(m)]]
        r=ok(m,pol,Ds); TL,E1,E2,E3,gap=terms(m,pol,Ds)
        allok&=r[0]
        print("m",m,"D*",Ds,"eps",r[1].str(5),"gap",r[2].str(5),"E1/TL",(E1/TL).str(3),"E2/TL",(E2/TL).str(3),"E3/TL",(E3/TL).str(3),"sqrt(D*) must exceed",mono(m,pol,Ds)[1].str(4),"ok",r[0],flush=True)
        print("   certificate at D*-1:",ok(m,pol,Ds-1)[0])
    print("tail certificate:","PASSED" if allok else "FAILED")
