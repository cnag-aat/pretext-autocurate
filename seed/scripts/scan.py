import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
def scan(rows_global, label, win=20, top=12, exclude=()):
    prof=np.asarray(M[rows_global,:],dtype=np.float32).mean(axis=0)
    N=PT.npix; oc=PT.orig_contig
    res=[]
    for c in range(PT.n_contigs):
        w=np.where(oc==c)[0]
        if len(w)==0: continue
        for s in range(w[0], w[-1]+1, win):
            e=min(s+win, w[-1]+1)
            if e-s<5: continue
            res.append((float(prof[s:e].mean()), PT.names[c], int(s-w[0]), int(e-w[0]), len(w)))
    res.sort(reverse=True)
    bg=np.median([r[0] for r in res])
    print('==',label,'  genome median window value: %.1f'%bg)
    for v,name,s,e,L in res[:top]:
        print('   %6.1f  %-16s [%d-%d] of %d px'%(v,name,s,e,L))
gA=lambda a,b: list(range(R['A'][0]+a,R['A'][0]+b))
gD=lambda a,b: list(range(R['D'][0]+a,R['D'][0]+b))
gC=lambda a,b: list(range(R['C'][0]+a,R['C'][0]+b))
gE=lambda a,b: list(range(R['E'][0]+a,R['E'][0]+b))
gB=lambda a,b: list(range(R['B'][0]+a,R['B'][0]+b))
gS=lambda a,b: list(range(R['S'][0]+a,R['S'][0]+b))
scan(gA(506,530)+gD(206,230), 'R1 start (A[506-530] + D[206-230])')
scan(gA(700,735)+gD(2,37), 'R3 end (A[700-735] + D[2-37])')
scan(gA(590,680)+gD(56,146), 'R3 middle')
scan(gA(493,503), 'A2 tail after rpt2 (A[493-503])')
scan(gE(322,350)+gS(14,42), 'X tail (E[322-350] + S[14-42])')
scan(gS(0,10), 'S[0-10]')
scan(gB(385,405), 'B end subtelomere B[385-405]')
scan(gC(1,20), 'C start subtelomere C[1-20]')
print('------ second batch')
scan(gA(450,480)+gC(340,371), 'Q end (A[450-480] + C[340-371])', top=14)
scan(gA(482,492), 'rpt2a A[482-492]', top=14)
scan(gC(372,381), 'C[372-381] (c, between Q and R2 in C)', top=10)
scan(gC(410,417), "C[410-417] (c', between rpt2 and R1 in C)", top=10)
scan(gA(100,112), 'A[100-112]', top=8)
scan(gA(118,148), 'A[118-148] (insertion vs B?)', top=10)
scan(gE(37,45), 'E[37-45] (E indel)', top=8)
scan(gB(36,47), 'B[36-47]', top=8)
