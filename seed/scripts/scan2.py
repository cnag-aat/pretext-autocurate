import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces2 import *
def scan(rows_global, label, win=20, top=12):
    prof=np.asarray(M[rows_global,:],dtype=np.float32).mean(axis=0)
    oc=PT.orig_contig; res=[]
    for c in range(PT.n_contigs):
        w=np.where(oc==c)[0]
        if len(w)==0: continue
        for s in range(w[0], w[-1]+1, win):
            e=min(s+win, w[-1]+1)
            if e-s<5: continue
            res.append((float(prof[s:e].mean()), PT.names[c], int(s-w[0]), int(e-w[0]), len(w)))
    res.sort(reverse=True)
    print('==',label,'  genome median: %.1f'%np.median([r[0] for r in res]))
    for v,name,s,e,L in res[:top]: print('   %6.1f  %-16s [%d-%d] of %d px'%(v,name,s,e,L))
gl=lambda k,a,b: list(range(R[k][0]+a,R[k][0]+b))
if __name__=='__main__':
    scan(gl('G',5,45),'G[5-45]')
    scan(gl('G',45,85),'G[45-85]')
    scan(gl('G',114,146),'G[114-146] (G-specific vs F)')
    scan(gl('F',290,312),'F[290-312]')
    scan(gl('H',125,165),'H[125-165]')
    scan(gl('H',170,205),'H[170-205]')
