import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces2 import *
import itertools
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def g(piece,pos): return R[piece][0]+pos
def plen(p): return R[p][1]-R[p][0]
def arm_contact(ax_, ay_, smin=1, smax=40, mat=M):
    n=smax-smin; tot=np.zeros((n,n)); 
    for (p1,j1,d1,L1) in ax_:
        for (p2,j2,d2,L2) in ay_:
            for a,s in enumerate(range(smin,smax)):
                x=j1+d1*s
                if s>=L1 or x<0 or x>=plen(p1): continue
                cs=[(b,j2+d2*t) for b,t in enumerate(range(smin,smax)) if t<L2 and 0<=j2+d2*t<plen(p2)]
                if not cs: continue
                tot[a,[b for b,_ in cs]]+=np.asarray(mat[g(p1,x),[g(p2,c) for _,c in cs]],dtype=float)
    return tot
def bands(T,smin):
    ii,jj=np.indices(T.shape); d=ii+jj+2*smin
    out=[]
    for lo,hi in ((0,20),(20,40),(40,60),(60,90)):
        m=(d>lo)&(d<=hi)
        out.append(T[m].mean() if m.any() else float('nan'))
    return out
arms={
 'Qe':[('A',481,-1,90),('C',372,-1,90)],
 'Ys':[('H',38,+1,200),('F',204,+1,85)],
 'R2R1':[('A',560,-1,56),('D',176,+1,56)],          # proximal of rpt2b: R2 then R1 (as in A3/D)
 'R3s':[('A',587,+1,150),('D',149,-1,150)],          # R3 distal of the first 20 px (not shared with H)
 'R3d':[('A',607,+1,130),('D',129,-1,130)],          # R3 distal of the F breakpoint (A606)
 'Gx':[('F',159,-1,150),('G',92,+1,20)],             # F continuation beyond its R3-like bit
}
keys=list(arms)
smin,smax=3,45
print('collapsed contact, by s+t band: <=20 | 20-40 | 40-60 | 60-90')
for x,y in itertools.combinations(keys,2):
    T=arm_contact(arms[x],arms[y],smin,smax)
    print('  %-10s'%(x+'·'+y),'  '.join('%7.1f'%v for v in bands(T,smin)))
