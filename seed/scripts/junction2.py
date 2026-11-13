import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
import itertools
def g(piece,pos): return R[piece][0]+pos
def plen(p): return R[p][1]-R[p][0]
def arm_contact(armsx, armsy, smin=3, smax=30, mat=M):
    """mean collapsed contact for s,t in [smin,smax): returns (T matrix summed over copy pairs)"""
    n=smax-smin
    tot=np.zeros((n,n)); cnt=np.zeros((n,n))
    for (p1,j1,d1) in armsx:
        for (p2,j2,d2) in armsy:
            for a,s in enumerate(range(smin,smax)):
                x=j1+d1*s
                if x<0 or x>=plen(p1): continue
                cs=[(b,j2+d2*t) for b,t in enumerate(range(smin,smax)) if 0<=j2+d2*t<plen(p2)]
                if not cs: continue
                bs=[b for b,_ in cs]; cols=[g(p2,c) for _,c in cs]
                tot[a,bs]+=np.asarray(mat[g(p1,x),cols],dtype=float)
    return tot
def bands(T,smin=3):
    ii,jj=np.indices(T.shape); d=ii+jj+2*smin
    return [T[(d>lo)&(d<=hi)].mean() for lo,hi in ((0,16),(16,30),(30,45),(45,60))]
def report(name, arms, smin=3, smax=30):
    print('==',name,' (bands of s+t: <=16, 16-30, 30-45, 45-60 px)')
    for x,y in itertools.combinations(arms.keys(),2):
        T=arm_contact(arms[x],arms[y],smin,smax)
        print('  %-6s'%(x+'·'+y), '  '.join('%7.1f'%v for v in bands(T,smin)))
J={'P':[('A',387,-1),('B',78,+1)],'Q':[('A',388,+1),('C',281,+1)],'U':[('B',77,-1),('E',164,-1)],'V':[('E',165,+1),('C',280,-1)]}
report('Junction J (A388 / B77 / C280 / E165)',J)
K={'W':[('E',126,-1),('C',136,-1)],'U':[('E',128,+1),('B',21,+1)],"V'":[('E',317,-1),('C',138,+1)],'X':[('E',319,+1),('B',13,-1),('S',11,+1)]}
report('Junction K (E127|E318 / C137 / B~17)',K)
