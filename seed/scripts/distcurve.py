import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces2 import *
def g(piece,pos): return R[piece][0]+pos
def plen(p): return R[p][1]-R[p][0]
def curve(ax_, ay_, offset, bins=((30,45),(45,60),(60,80),(80,100),(100,130),(130,170)), smax=120, occ=False):
    """mean collapsed contact between arm x at s and arm y at t, binned by real distance s+t+offset.
    collapsed = sum over copy pairs; normalised per (s,t) cell."""
    acc={b:[0.0,0] for b in bins}
    S=range(1,smax)
    # build collapsed matrix
    T=np.zeros((smax,smax)); V=np.zeros((smax,smax),bool)
    for (p1,j1,d1,L1) in ax_:
        for (p2,j2,d2,L2) in ay_:
            for s in S:
                x=j1+d1*s
                if s>=L1 or x<0 or x>=plen(p1): continue
                ts=[t for t in S if t<L2 and 0<=j2+d2*t<plen(p2)]
                if not ts: continue
                vals=np.asarray(M[g(p1,x),[g(p2,j2+d2*t) for t in ts]],dtype=float)
                if occ: vals=(vals>0).astype(float)
                T[s,ts]+=vals; V[s,ts]=True
    ii,jj=np.indices(T.shape); d=ii+jj+offset
    out=[]
    for lo,hi in bins:
        m=V&(d>=lo)&(d<hi)
        out.append(T[m].mean() if m.sum()>20 else float('nan'))
    return out
arms={
 'Qe':[('A',481,-1,93),('C',372,-1,91)],
 'Ys':[('H',38,+1,84),('F',204,+1,84)],
 'Gx':[('F',159,-1,150)],
 'P':[('A',387,-1,300),('B',78,+1,300)],'Q':[('A',388,+1,92),('C',281,+1,90)],'U':[('B',77,-1,55),('E',164,-1,37)],'V':[('E',165,+1,150),('C',280,-1,140)],
 'W':[('E',126,-1,120),('C',136,-1,110)],"V'":[('E',317,-1,150),('C',138,+1,140)],
 'refL':[('A',249,-1,140),('B',194,+1,80)],'refR':[('A',250,+1,80),('B',193,-1,80)],
 'refL2':[('E',239,-1,70),('C',206,+1,70)],'refR2':[('E',240,+1,75),('C',205,-1,65)],
 'R3a':[('A',650,-1,80),('D',86,+1,80)],'R3b':[('A',651,+1,85),('D',85,-1,85)],
}
print('distance bins (px):            30-45   45-60   60-80  80-100 100-130 130-170')
def show(lab,x,y,off):
    c=curve(arms[x],arms[y],off); print('%-28s'%lab,'  '.join('%6.1f'%v for v in c))
show('ref cis P interior (2 copies)','refL','refR',0)
show('ref cis V interior (2 copies)','refL2','refR2',0)
show('ref cis R3 interior (A+D)','R3a','R3b',0)
show('J: P·Q','P','Q',0); show('J: P·U','P','U',0); show('J: Q·V','Q','V',0); show('J: U·V','U','V',0)
show('J: P·V (neg)','P','V',0); show('J: Q·U (neg)','Q','U',0)
show("K: W·V'",'W',"V'",0)
show('L: Q·Y  (+22)','Qe','Ys',22)
show('L: G·Y  (+45)','Gx','Ys',45)
show('L: Q·G (neg)','Qe','Gx',0)
