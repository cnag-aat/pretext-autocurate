import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lo,hi=380,737
apos=np.arange(lo,hi)
def interp_track(piece, anchors):
    """anchors: list of segments [(a0,p0,a1,p1)] mapping A coord -> piece coord linearly; returns global idx array (or -1)"""
    out=-np.ones(len(apos),dtype=int)
    for a0,p0,a1,p1 in anchors:
        for a in range(a0,a1+1):
            if lo<=a<hi:
                p=int(round(p0+(p1-p0)*(a-a0)/max(1,(a1-a0))))
                out[a-lo]=R[piece][0]+p
    return out
tA=np.array([R['A'][0]+a for a in apos])
tC=interp_track('C',[(386,280,421,310),(422,314,472,361),(473,366,479,371),(502,458,545,418),(546,384,568,402)])
tD=interp_track('D',[(504,232,736,0)])
tracks={'A':tA,'C':tC,'D':tD}
n=len(apos)
def collapsed(pairs=None, mat=M):
    tot=np.zeros((n,n))
    for k1,t1 in tracks.items():
        for k2,t2 in tracks.items():
            if pairs is not None and (k1,k2) not in pairs: continue
            v1=np.where(t1>=0)[0]; v2=np.where(t2>=0)[0]
            sub=np.asarray(mat[np.ix_(t1[v1],t2[v2])],dtype=float)
            tot[np.ix_(v1,v2)]+=sub
    return tot
T=collapsed()
np.save('/home/claude/work/data/collapsedR.npy',T)
fig,axs=plt.subplots(1,2,figsize=(22,11))
for ax,(TT,title) in zip(axs,[(T,'collapsed over all copies (A + C + D), A coordinates'),(collapsed({('A','A')}),'A x A only')]):
    ax.imshow((np.clip(TT,0,500)/500)**0.45,cmap='magma_r',vmin=0,vmax=1,extent=(lo,hi,hi,lo),interpolation='nearest')
    for x,c in [(388,'g'),(480,'b'),(503.5,'purple'),(545,'c'),(569,'c')]:
        ax.axvline(x,color=c,lw=0.5); ax.axhline(x,color=c,lw=0.5)
    ax.set_title(title); ax.set_xticks(range(380,737,20)); ax.set_yticks(range(380,737,20)); ax.tick_params(labelsize=7)
fig.savefig('/home/claude/work/out/collapsedR.png',dpi=80,bbox_inches='tight')
# virtual 4C profiles
def v4c(rows, label):
    prof=T[[r-lo for r in rows],:].mean(axis=0)
    return prof
for label,rows in [('Q end A[455-478]',range(455,479)),('R1 A[508-540]',range(508,541)),('R2 A[547-560]',range(547,561)),('R3 start A[572-600]',range(572,601)),('R3 mid A[640-680]',range(640,681))]:
    prof=v4c(rows,label)
    bins=[(390,420),(420,450),(450,480),(482,492),(493,503),(504,515),(515,530),(530,545),(545,561),(561,569),(570,590),(590,620),(620,660),(660,700),(700,737)]
    print('%-22s'%label, ' '.join('%s:%5.1f'%('%d-%d'%b, prof[b[0]-lo:b[1]-lo].mean()) for b in bins))
