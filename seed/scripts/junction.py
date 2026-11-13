import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
# arm copies around junction J: list of (piece, J_coordinate, direction) ; position at distance s from J = J + dir*s
arms={
 'P':[('A',387,-1),('B',78,+1)],
 'Q':[('A',388,+1),('C',281,+1)],
 'U':[('B',77,-1),('E',164,-1)],
 'V':[('E',165,+1),('C',280,-1)],
}
def g(piece,pos): return R[piece][0]+pos
def arm_contact(x,y,smin=4,smax=40,mat=M):
    """collapsed contact matrix between arm x (distance s) and arm y (distance t) summed over all copy pairs"""
    n=smax-smin
    tot=np.zeros((n,n))
    for (p1,j1,d1) in arms[x]:
        for (p2,j2,d2) in arms[y]:
            r=[g(p1,j1+d1*s) for s in range(smin,smax)]
            c=[g(p2,j2+d2*t) for t in range(smin,smax)]
            tot+=np.asarray(mat[np.ix_(r,c)],dtype=float)
    return tot
import itertools
print('MQ0 collapsed contact between arms near junction J (mean intensity, by distance band s+t)')
print('pair   s+t<=20   20-40   40-60   60-80')
res={}
for x,y in itertools.combinations('PQUV',2):
    T=arm_contact(x,y)
    n=T.shape[0]; ii,jj=np.indices(T.shape); d=ii+jj+8
    vals=[T[(d>lo)&(d<=hi)].mean() for lo,hi in ((0,20),(20,40),(40,60),(60,80))]
    res[x+y]=vals
    print(x+y, '  ', '  '.join('%7.1f'%v for v in vals))
# reference: true cis contact at same distances inside clean single-junction-free regions, collapsed over 2 copies
# use P interior: A[200..] vs B copy: pseudo junction at A=250 (B ~ 193)
def ref(pA,jA,pB,jB,dB):
    arms['L']=[(pA,jA-1,-1),(pB,jB-dB,-dB)]
    arms['Rr']=[(pA,jA,+1),(pB,jB,dB)]
    T=arm_contact('L','Rr'); ii,jj=np.indices(T.shape); d=ii+jj+8
    return [T[(d>lo)&(d<=hi)].mean() for lo,hi in ((0,20),(20,40),(40,60),(60,80))]
print('ref cis (P interior, A250/B193):', ['%.1f'%v for v in ref('A',250,'B',193,-1)])
print('ref cis (P interior, A200/B240):', ['%.1f'%v for v in ref('A',200,'B',240,-1)])
print('ref cis (V interior, E240/C205):', ['%.1f'%v for v in ref('E',240,'C',205,-1)])
print('ref cis (W interior, E80/C92):', ['%.1f'%v for v in ref('E',80,'C',92,+1)])
# background: unrelated arms
arms['Z']=[('D',100,+1),('D',99,-1)]
for x in 'PQUV':
    T=arm_contact(x,'Z'); print('bg',x,'x D: %.1f'%T.mean())
