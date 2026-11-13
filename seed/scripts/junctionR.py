import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
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
arms={
 'Qe':[('A',480,-1,90),('C',372,-1,90)],
 'R1s':[('A',503,+1,40),('D',233,-1,40),('C',459,-1,40)],
 'R1e':[('A',545,-1,40),('D',191,+1,40),('C',417,+1,40)],
 'R2s':[('A',545,+1,16),('D',191,-1,16),('C',382,+1,18)],
 'R2e':[('A',561,-1,16),('D',176,+1,16),('C',401,-1,18)],
 'R3s':[('A',569,+1,160),('D',165,-1,160)],
}
keys=list(arms)
print('collapsed mean contact, s,t in 2..15 px from each arm end')
print('        '+'  '.join('%6s'%k for k in keys))
for x in keys:
    row=[]
    for y in keys:
        if x==y: row.append('   -  '); continue
        T=arm_contact(arms[x],arms[y],2,15)
        row.append('%6.1f'%T.mean())
    print('%6s  '%x+'  '.join(row))
n=len(keys)
fig,axs=plt.subplots(n,n,figsize=(2.6*n,2.6*n))
for i,x in enumerate(keys):
    for j,y in enumerate(keys):
        ax=axs[i,j]
        if i==j: ax.axis('off'); continue
        T=arm_contact(arms[x],arms[y],1,40)
        ax.imshow((np.clip(T,0,400)/400)**0.5,cmap='magma_r',vmin=0,vmax=1,origin='lower',extent=(1,40,1,40))
        ax.set_title(f'{x} (y) x {y} (x)',fontsize=8); ax.tick_params(labelsize=5)
fig.savefig('/home/claude/work/out/cornerR.png',dpi=85,bbox_inches='tight')
