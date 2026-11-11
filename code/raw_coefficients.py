import os
HERE=os.path.dirname(os.path.abspath(__file__))
import json
for m in (1,2,3,5,7):
    d={tuple(map(int,k.split(','))):v for k,v in json.load(open(os.path.join(HERE,'data',f'cF_{m}.json'))).items()}
    F={}
    for (n,l),v in d.items():
        D=4*m*n-l*l; F.setdefault(D,set()).add((-1)**(l+1)*v)
    assert all(len(s)==1 for s in F.values()); F={D:s.pop() for D,s in F.items()}
    def G(D):
        l=[l for l in range(m+1) if (D+l*l)%(4*m)==0][0]
        if (D+l*l)//(4*m)<-1: return 0
        return F[D]
    bad=[]
    for J in range(0,300):
        for n in range(m,400):
            D=4*m*n-J*J
            if D<=0 or D>1000: continue
            Dp=D-4*(J+1)
            if G(D)-G(Dp)<=0: bad.append((n,J,D,Dp,G(D),G(Dp)))
    print(m,'polar F',sorted((D,v) for D,v in F.items() if D<=0),'bad',len(bad),bad[:6])
