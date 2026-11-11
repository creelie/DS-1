"""Write the tables of Appendix A (phi_m and the polar terms of phi_m/Delta) as LaTeX.
Usage: python3 make_appendix.py > appendix_tables.tex   (a few minutes)."""
import os, json
from fractions import Fraction as Fr
HERE=os.path.dirname(os.path.abspath(__file__))
import identity
def frac(x):
    x=Fr(str(x))
    if x.denominator==1: return f"{x.numerator}"
    s='-' if x<0 else ''
    return f"{s}\\tfrac{{{abs(x.numerator)}}}{{{x.denominator}}}"
def mono(a,b,i,j):
    t=[]
    for name,e in (("E_4",a),("E_6",b),("A",i),("B",j)):
        if e==1: t.append(name)
        elif e>1: t.append(f"{name}^{{{e}}}")
    return " ".join(t)
def phi_table(m):
    cF=identity.fast.compute(m,16)
    mons,x,rank,ncols,ok,nk=identity.decompose(m,14,cF)
    assert ok and rank==ncols
    rows=[(mono(*mo),frac(xi)) for mo,xi in zip(mons,x[:-1])]
    h=(len(rows)+1)//2
    left,right=rows[:h],rows[h:]
    out=[f"\\noindent$m={m}$:\\nopagebreak\\\\[2pt]","{\\small\\begin{tabular}{@{}l@{\\ \\ }l@{\\qquad}l@{\\ \\ }l@{}}"]
    for k in range(h):
        L=left[k]; R=right[k] if k<len(right) else ("","")
        out.append(f"${L[0]}$ & ${L[1]}$ & " + (f"${R[0]}$ & ${R[1]}$" if R[0] else "&") + " \\\\")
    out.append("\\end{tabular}}\\par\\medskip")
    return "\n".join(out)
def polar_table(m,pol):
    rows=[f"$({D},{l})$ & ${a}$" for D,l,a in pol]
    ncol=3 if m<=13 else 2
    size="\\small" if m<=13 else "\\scriptsize"
    head="{"+size+"\\begin{tabular}{@{}"+"@{\\qquad}".join(["r@{\\ \\ }r"]*ncol)+"@{}}"
    out=[f"\\noindent$m={m}$:\\nopagebreak\\\\[2pt]"]
    # long tables are set in blocks of at most 30 lines, read column by column, so that they can break across pages
    for b in range(0,len(rows),30*ncol):
        blk=rows[b:b+30*ncol]; h=-(-len(blk)//ncol)
        cols=[blk[i*h:(i+1)*h] for i in range(ncol)]
        out.append(head)
        for k in range(h):
            cells=[cols[i][k] if k<len(cols[i]) else "&" for i in range(ncol)]
            out.append(" & ".join(cells)+" \\\\")
        out.append("\\end{tabular}}\\par")
    out[-1]+="\\medskip"
    return "\n".join(out)
if __name__=="__main__":
    pols=json.load(open(os.path.join(HERE,'data','polar.json')))
    print("%% phi_m tables")
    for m in (3,5,7,11,13): print(phi_table(m))
    print("%% polar tables")
    from classnum import MSET
    for m in MSET: print(polar_table(m,pols[str(m)]))
