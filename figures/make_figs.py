"""Figures of the paper.  Run from the repository root: python3 figures/make_figs.py
Writes PNG files into paper/."""
import os, sys, json, glob, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(ROOT,'code'))
OUT=os.path.join(ROOT,'paper')
plt.rcParams.update({"font.size":10,"font.family":"serif","mathtext.fontset":"cm"})
cF={}
for f in glob.glob(os.path.join(ROOT,'code','data','cF_*.json')):
    m=int(f.split('_')[-1][:-5]); cF[m]={tuple(map(int,k.split(','))):v for k,v in json.load(open(f)).items()}
def reduce(a,b,c):
    while True:
        if c>a: a,c=c,a
        if abs(b)>c:
            k=(b+c)//(2*c); a,b=a-b*k+c*k*k,b-2*c*k; continue
        return a,abs(b),c
def dstar(n,l,m):
    D=4*m*n-l*l
    if D<=0: return 0
    a,b,c=reduce(n,l,m); return (-1)**(b+1)*cF[c][(a,b)]
def Fm(m):
    F={}
    for (n,l),v in cF[m].items():
        D=4*m*n-l*l
        if 0<=l<=m and D>0: F[D]=(-1)**(l+1)*v
    return F
MS=[1,2,3,5,7,11,13]
cols=plt.cm.viridis(np.linspace(0,0.9,len(MS)))

# 1. growth of the single-centred index F_m(D) and the leading Bessel term
from flint import arb
fig,ax=plt.subplots(figsize=(6.2,3.8))
for m,c in zip(MS,cols):
    F=Fm(m)
    red={4*m*n-l*l for (n,l) in cF[m] if 0<=l<=m and n>=m}   # discriminants of reduced charges
    Ds=sorted(D for D in F if D<=1000 and D in red)
    ax.plot(Ds,[math.log10(float(F[D])) for D in Ds],'.',ms=2.2,color=c,label=f"$m={m}$")
    al=math.pi*math.sqrt((m+4)/m); nu=11.5
    xs=np.linspace(20,1000,300)
    T=[math.log10(2*math.pi/math.sqrt(2*m)*(m+2))+ (nu/2)*math.log10(m*(m+4)/x)+ float((arb(al)*arb(x).sqrt()).bessel_i(arb(nu)).log().mid())/math.log(10) for x in xs]
    ax.plot(xs,T,'-',lw=0.7,color=c)
ax.set_xlabel("discriminant $D=4mn-l^2$"); ax.set_ylabel(r"$\log_{10}F_m(D)$")
ax.legend(ncol=4,fontsize=8,frameon=False,loc="upper left")
ax.set_title("Single-centred index (dots) and leading Bessel term (lines)",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_growth.png"),dpi=200); plt.close(fig)

# 2. tail certificate: eps(D) against gap(D)/2
import tail
pols=json.load(open(os.path.join(ROOT,'code','data','polar.json')))
DST={1:680,2:371,3:292,5:258,7:313,11:720,13:1030}
fig,axs=plt.subplots(1,2,figsize=(7.2,3.2))
for m,c in zip(MS,cols):
    pol=[tuple(t) for t in pols[str(m)]]
    Ds=np.arange(max(60,DST[m]//3),2*DST[m]+1,max(4,DST[m]//60))
    e=[];g=[]
    for D in Ds:
        TL,E1,E2,E3,gap=tail.terms(m,pol,int(D),DST[m])
        e.append(float(((E1+E2+E3)/TL).mid())); g.append(float(gap.mid())/2)
    axs[0].semilogy(Ds,e,'-',color=c,lw=1,label=f"$m={m}$"); axs[0].semilogy(Ds,g,'--',color=c,lw=0.8)
    axs[0].plot([DST[m]],[float(tail.ok(m,pol,DST[m])[1].mid())],'o',color=c,ms=3)
    axs[1].plot(Ds/DST[m],np.array(e)/np.array(g),'-',color=c,lw=1)
axs[0].set_xlabel("$D$"); axs[0].set_ylabel(r"$\varepsilon_m(D)$ (solid), $\gamma_m(D)/2$ (dashed)")
axs[0].set_ylim(1e-2,1e2); axs[0].legend(fontsize=7,frameon=False,ncol=2)
axs[1].axhline(1,color='k',lw=0.6); axs[1].axvline(1,color='k',lw=0.6,ls=':')
axs[1].set_xlabel(r"$D/D^*_m$"); axs[1].set_ylabel(r"$2\varepsilon_m(D)/\gamma_m(D)$"); axs[1].set_ylim(0,3)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_tail.png"),dpi=200); plt.close(fig)

# 3. relative margins of the fixed-J index over all charges with D <= 1000
pts=[]
for J in range(0,252):
    for D in range(1,1001):
        if (D+J*J)%4: continue
        P=(D+J*J)//4
        for m in range(1,P+1):
            if P%m: continue
            n=P//m
            d1=dstar(n,J,m); d2=dstar(n,J+2,m)
            pts.append((D,(d1-d2)/d1,m))
pts=np.array(pts)
fig,ax=plt.subplots(figsize=(6.2,3.4))
sc=ax.scatter(pts[:,0],pts[:,1],c=np.log(pts[:,2]),s=0.6,cmap="plasma",rasterized=True)
ax.set_xlabel("$D$"); ax.set_ylabel(r"$d(n,J,m)/d^*(n,J,m)$")
ax.axhline(0,color='k',lw=0.5); ax.set_ylim(-0.05,1.05)
cb=fig.colorbar(sc,ax=ax); cb.set_label(r"$\log m$")
ax.set_title(f"{len(pts)} charges, minimum ratio {pts[:,1].min():.4f}",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_margins.png"),dpi=200); plt.close(fig)

# 4. raw coefficients against the reduced prescription (m=5)
m=5; F=Fm(m)
fig,ax=plt.subplots(figsize=(6.2,3.2))
rawp=[];redp=[]
for n in range(m,40):
    for J in range(0,60):
        D=4*m*n-J*J; Dp=D-4*(J+1)
        if D<=0 or D>1000: continue
        def G(DD):
            l=[l for l in range(m+1) if (DD+l*l)%(4*m)==0][0]
            return 0 if (DD+l*l)//(4*m)<-1 else F.get(DD,0)
        raw=F[D]-G(Dp); red=dstar(n,J,m)-dstar(n,J+2,m)
        rawp.append((D,np.sign(raw)*math.log10(1+abs(raw)))); redp.append((D,math.log10(1+red)))
rawp=np.array(rawp); redp=np.array(redp)
ax.plot(redp[:,0],redp[:,1],'.',ms=2,color='tab:blue',label="index $d(n,J,5)$ (reduced charges)")
neg=rawp[rawp[:,1]<0]
ax.plot(neg[:,0],neg[:,1],'x',ms=4,color='tab:red',label="raw differences (no reduction), negative values")
ax.axhline(0,color='k',lw=0.5)
ax.set_xlabel("$D$"); ax.set_ylabel(r"$\mathrm{sign}\cdot\log_{10}(1+|\cdot|)$"); ax.legend(fontsize=8,frameon=False)
ax.set_title("$m=5$, $n\\geq m$: raw differences versus the index $d(n,J,5)$",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_raw.png"),dpi=200); plt.close(fig)

# 5. 3D view of the fixed-J index for m = 7 over (n, J)
m=7
fig=plt.figure(figsize=(6.0,4.4)); ax=fig.add_subplot(111,projection='3d')
X=[];Y=[];Z=[]
for n in range(1,37):
    for J in range(0,2*int(math.sqrt(m*n))+1):
        D=4*m*n-J*J
        if D<=0 or D>1000: continue
        d=dstar(n,J,m)-dstar(n,J+2,m)
        X.append(n);Y.append(J);Z.append(math.log10(d))
ax.scatter(X,Y,Z,c=Z,cmap='viridis',s=6,depthshade=True)
ax.set_xlabel("$n$"); ax.set_ylabel("$J$"); ax.set_zlabel(r"$\log_{10} d(n,J,7)$")
ax.view_init(elev=24,azim=-58)
ax.set_title("Fixed angular momentum index, $Q_1Q_5=7$",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_index3d.png"),dpi=200); plt.close(fig)

# 6. poles: chi/A has double poles at the lattice points; the Appell-Lerch sum has the same principal parts
fig,ax=plt.subplots(figsize=(5.2,3.6))
tau=complex(0.32,0.9)
for s in range(-2,3):
    for mu in range(-2,3):
        z=s*tau+mu
        ax.plot(z.real,z.imag,'o',color='tab:red' if s==0 else 'tab:purple',ms=5)
for s in range(-2,3):
    ax.plot([s*tau.real-3,s*tau.real+3],[s*tau.imag]*2,color='0.85',lw=0.6,zorder=0)
ax.fill_between([-3,3],[0]*2,[tau.imag]*2,color='tab:blue',alpha=0.08)
ax.text(1.55,0.42,r"$|q|<|\zeta|<1$",fontsize=9)
bb=dict(boxstyle='round,pad=0.2',fc='white',ec='0.7',lw=0.5)
ax.annotate(r"$z=0$: $\chi_{m+1}/A\sim p_{24}(m+1)/(2\pi i z)^2$",xy=(0,0),xytext=(-2.9,1.35),fontsize=8,bbox=bb,arrowprops=dict(arrowstyle='->',lw=0.6))
ax.annotate(r"$z=-\tau$: principal part at $0$ times $q^{-m}e(2m\varepsilon)$",xy=(-tau.real,-tau.imag),xytext=(-2.9,-2.25),fontsize=8,bbox=bb,arrowprops=dict(arrowstyle='->',lw=0.6))
ax.set_xlim(-3,3); ax.set_ylim(-2.4,2.2); ax.set_aspect('equal')
ax.set_xlabel(r"$\mathrm{Re}\,z$"); ax.set_ylabel(r"$\mathrm{Im}\,z$")
ax.set_title(r"Double poles of $\psi_m$ in the $z$-plane",fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_poles.png"),dpi=200); plt.close(fig)
print("figures written to",OUT)
