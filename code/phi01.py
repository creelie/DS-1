# c(D) of 2 phi_{0,1}, via phi_{0,1}=4 sum_i (th_i(z)/th_i(0))^2 ; only y^0 and y^1 coefficients needed
import sys,json,os
from flint import fmpq_poly, fmpq
NQ=int(sys.argv[1]); E=8*NQ+8   # exponents in q^{1/8}
def ser(f):  # f: dict exp->coef  -> fmpq_poly truncated
    a=[0]*(E+1)
    for e,v in f.items():
        if e<=E: a[e]+=v
    return fmpq_poly(a)
def trunc(p): 
    c=p.coeffs()[:E+1]; return fmpq_poly(c)
R=int((E)**0.5)+3
res={}
for name in ('2','3','4'):
    # coefficient of y^L in th(z)^2 for L=0,1 ; and th(0)^2
    if name=='3': terms=[(n, 4*n*n, 1) for n in range(-R,R+1)]           # (y-exp in 1/2 units*... use 2n half-units)
    if name=='4': terms=[(n, 4*n*n, (-1)**(n%2)) for n in range(-R,R+1)]
    if name=='2': terms=[(n, (2*n+1)**2, 1) for n in range(-R,R)]
    # y-exponent in half units: th3/th4: 2n ; th2: 2n+1
    def yh(n): return 2*n if name!='2' else 2*n+1
    sq={0:{},2:{}}  # half-unit y exponent target 0 (y^0) and 2 (y^1)
    zero={}
    for (a,ea,sa) in terms:
        for (b,eb,sb) in terms:
            e=ea+eb
            if e>E: continue
            zero[e]=zero.get(e,0)+sa*sb
            t=yh(a)+yh(b)
            if t in sq: sq[t][e]=sq[t][e]+sa*sb if e in sq[t] else sa*sb
    z=ser(zero); m=min(k for k,v in zero.items() if v)
    # 1/z : z = q^m * u
    u=fmpq_poly(z.coeffs()[m:]); 
    def tr(p,k): return fmpq_poly(p.coeffs()[:k])
    inv=fmpq_poly([fmpq(1)/u.coeffs()[0]]); k=1
    while k<E+1:
        k=min(2*k,E+1)
        inv=tr(inv*(2-tr(u*inv,k)),k)
    for t in (0,2):
        num=ser(sq[t]); 
        numc=num.coeffs()
        # divide by q^m
        assert all(c==0 for c in numc[:m])
        nn=fmpq_poly(numc[m:])
        p=trunc(nn*inv)
        res.setdefault(t,fmpq_poly(0)); res[t]=res[t]+8*p
c={}
for t,l in ((0,0),(2,1)):
    co=res[t].coeffs()
    for e,v in enumerate(co):
        if v==0: continue
        assert e%8==0 and v.q==1,(e,v)
        n=e//8
        if n<=NQ-2: c[4*n-l*l]=int(v.p)
json.dump({str(k):v for k,v in sorted(c.items())},open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','c2phi01.json'),'w'))
print(len(c), [ (k,c[k]) for k in sorted(c)[:12]], max(c))
