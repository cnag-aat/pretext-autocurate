"""Compare a tool savestate with a curator savestate of the same map.

Usage: python3 tools/compare_layouts.py MAP.pretext TOOL.savestate CURATOR.savestate
Prints junction agreement both ways, curator chromosome -> tool paint groups, and how many
curator-tagged pieces sit in the right group. Exploratory script, not part of the package.
"""
import sys,glob,collections
import os
R=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/'
sys.path.insert(0,R+'src')
import numpy as np
from pretext_autocurate import *
M=sys.argv[1]
pt=Pretext(M,load_graphs=False); N=pt.npix; bpp=pt.bp_per_pix
def load(path):
    st=SaveState.load(path); ms=MapState(pt).replay(st.edits); paint=dict(st.scaffs); fl=dict(st.meta_flags)
    pos=np.empty(N,int); pos[ms.pix]=np.arange(N); ch=np.zeros(N,int); tg=np.zeros(N,np.int64)
    for i,c in enumerate(ms.contigs):
        px=ms.pix[c.start:c.start+c.length]; ch[px]=paint.get(i,0); tg[px]=fl.get(i,0)&3
    return st,ms,pos,ch,tg
tool=load(sys.argv[2]); cur=load(sys.argv[3])
def junctions(L,minpx=15):
    st,ms,pos,ch,tg=L; out=[]
    big=[c for c in ms.contigs if c.length>=minpx and ch[ms.pix[c.start]]>0 and tg[ms.pix[c.start]]==0]
    for a,b in zip(big,big[1:]):
        k=min(6,a.length//2,b.length//2)      # compare a few pixels inside each piece: boundaries differ by a pixel or two
        pa=int(ms.pix[a.start+a.length-1-k]); pb=int(ms.pix[b.start+k])
        if ch[pa]==ch[pb]: out.append((pa,pb))
    return out
def supported(j,L,D=80):
    st,ms,pos,ch,tg=L; pa,pb=j
    return ch[pa]==ch[pb] and ch[pa]>0 and abs(pos[pa]-pos[pb])<=D
loc=lambda p:'%s:%.2f'%(pt.names[pt.orig_scaffold[p]],pt.rel_coord[p]*bpp/1e6)
jt=junctions(tool); jc=junctions(cur)
orig=lambda j: abs(j[0]-j[1])<=13
print('tool junctions between untagged contigs>=15px:',len(jt),'(new, not original adjacency: %d)'%sum(not orig(j) for j in jt))
bad=[j for j in jt if not supported(j,cur)]
print('  supported by curation: %d / %d'%(len(jt)-len(bad),len(jt)))
for j in bad: print('     not in curation:',loc(j[0]),'|',loc(j[1]))
miss=[j for j in jc if not supported(j,tool)]
print('curator junctions:',len(jc),'(new: %d)'%sum(not orig(j) for j in jc),' reproduced by tool: %d'%(len(jc)-len(miss)))
for j in miss: print('     missing:',loc(j[0]),'|',loc(j[1]))
# chromosome purity: for each curator chromosome, number of tool paint ids covering its untagged bp
st,ms,pos,ch,tg=cur; tch=tool[3]
print('curator chromosome -> tool groups (untagged sequence, groups with >=1% of the chromosome):')
for c in sorted(set(ch[ch>0])):
    m=(ch==c)&(tg==0)
    if not m.any(): continue
    cnt=collections.Counter(tch[m].tolist()); tot=m.sum()
    parts=[(k,v) for k,v in cnt.most_common() if v>=0.01*tot]
    print('  chr %2d %6.1f Mb: %s'%(c,tot*bpp/1e6,', '.join('%s %.0f%%'%('unpainted' if k==0 else 'g%d'%k,100*v/tot) for k,v in parts)))
# end attachment: tagged pieces on the same tool group as the main part of their curator chromosome
main={c:collections.Counter(tch[(ch==c)&(tg==0)].tolist()).most_common(1)[0][0] for c in set(ch[ch>0]) if ((ch==c)&(tg==0)).any()}
ok=bad=0
for i,cc in enumerate(ms.contigs):
    p=int(ms.pix[cc.start])
    if tg[p] and ch[p]>0:
        if ch[p] not in main: continue
        if tch[p]==main[ch[p]]: ok+=1
        else: bad+=1
print('curator unloc/haplotig pieces painted into the tool group holding most of their chromosome: %d, elsewhere/unpainted: %d'%(ok,bad))
