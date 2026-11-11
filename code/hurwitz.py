from fractions import Fraction as Fr
def Hcl(N):
    if N==0: return Fr(-1,12)
    if N%4 in (1,2): return Fr(0)
    s=Fr(0); b=N%2
    while 3*b*b<=N:
        if (b*b+N)%4==0:
            ac=(b*b+N)//4
            a=max(b,1)
            while a*a<=ac:
                if ac%a==0:
                    c=ac//a
                    if a==b==c: s+=Fr(1,3)
                    elif b==0 and a==c: s+=Fr(1,2)
                    elif b==0 or a==b or a==c: s+=1
                    else: s+=2
                a+=1
        b+=2
    return s
