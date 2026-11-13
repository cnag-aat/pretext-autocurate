import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
from scipy.ndimage import uniform_filter1d, convolve

def diag_scores(X, w=4):
    """forward (i+k,j+k) and reverse (i+k,j-k) diagonal means over window 2w+1."""
    h,wd=X.shape
    F=np.zeros_like(X); Rv=np.zeros_like(X); n=2*w+1
    for k in range(-w,w+1):
        # forward
        src=np.zeros_like(X)
        i0,i1=max(0,-k),min(h,h-k); j0,j1=max(0,-k),min(wd,wd-k)
        src[i0:i1,j0:j1]=X[i0+k:i1+k,j0+k:j1+k]; F+=src
        src=np.zeros_like(X)
        j0,j1=max(0,k),min(wd,wd+k)
        src[i0:i1,j0:j1]=X[i0+k:i1+k,j0-k:j1-k]; Rv+=src
    return F/n, Rv/n

def trace(p,q,w=4,thr=70,minlen=6, selfmask=30):
    X=block(p,q)
    if p==q:
        h=X.shape[0]
        ii,jj=np.indices(X.shape); X=np.where(np.abs(ii-jj)<selfmask,0,X)
    F,Rv=diag_scores(X,w)
    segs=[]
    for name,S,sgn in (('+',F,1),('-',Rv,-1)):
        best=S.argmax(axis=1); val=S.max(axis=1)
        # chain consecutive rows where best col moves by sgn (+-2)
        cur=None
        for i in range(S.shape[0]):
            ok = val[i]>=thr
            if ok and cur is not None and abs((best[i]-cur[-1][1]) - sgn)<=2 and i-cur[-1][0]<=2:
                cur.append((i,int(best[i]),float(val[i])))
            else:
                if cur is not None and len(cur)>=minlen: segs.append((name,cur))
                cur=[(i,int(best[i]),float(val[i]))] if ok else None
        if cur is not None and len(cur)>=minlen: segs.append((name,cur))
    out=[]
    for name,c in segs:
        out.append((p,q,name,c[0][0],c[-1][0],c[0][1],c[-1][1],len(c),round(float(np.mean([v for _,_,v in c])),0)))
    return out

if __name__=='__main__':
    allsegs=[]
    ks='ABCDES'
    for i,p in enumerate(ks):
        for q in ks[i:]:
            for s in trace(p,q):
                allsegs.append(s)
                print('%s[%d-%d] ~ %s[%d-%d] (%s) n=%d score=%d'%(s[0],s[3],s[4],s[1],s[5],s[6],s[2],s[7],s[8]))
