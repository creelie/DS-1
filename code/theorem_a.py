import os
HERE=os.path.dirname(os.path.abspath(__file__))
import json,sys,glob
Dmax=int(sys.argv[1])
cF={}
for f in glob.glob(os.path.join(HERE,'data','cF_*.json')):
    m=int(f.split('_')[-1][:-5]); cF[m]={tuple(map(int,k.split(','))):v for k,v in json.load(open(f)).items()}
def reduce(a,b,c):
    """Gauss reduction of the positive definite form [a,b,c] = a x^2 + b xy + c y^2.
    Returns (a',|b'|,c') with |b'| <= c' <= a' in the same SL2(Z) class."""
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
    tab=cF.get(c)
    if tab is None or (a,b) not in tab: raise KeyError((n,l,m,'->',a,b,c))
    return (-1)**(b+1)*tab[(a,b)]
# 4D: all reduced forms with 0<D<=Dmax
bad4=[];n4=0
for c in range(1,int((Dmax/3)**0.5)+2):
    for b in range(0,c+1):
        for a in range(c, (Dmax+b*b)//(4*c)+1):
            D=4*a*c-b*b
            if D<=0 or D>Dmax: continue
            n4+=1
            v=dstar(a,b,c)
            if v<=0: bad4.append((a,b,c,D,v))
print("reduced 4D forms",n4,"nonpositive",bad4[:10])
# 5D
bad5=[];n5=0;worst=(9,None)
for J in range(0,Dmax//4+2):
    for D in range(1,Dmax+1):
        if (D+J*J)%4: continue
        P=(D+J*J)//4
        for m in range(1,P+1):
            if P%m: continue
            n=P//m
            d1=dstar(n,J,m); d2=dstar(n,J+2,m)
            n5+=1
            if d1-d2<=0: bad5.append((m,n,J,D,d1,d2))
            elif d2>0:
                r=(d1-d2)/d1
                if r<worst[0]: worst=(r,(m,n,J,D))
print("5D charges (m,n,J)",n5,"nonpositive",bad5[:10],"min relative margin",worst)
