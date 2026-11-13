import sys, struct, os; sys.path.insert(0,'/home/claude/work/scripts')
from pretext import *
import numpy as np
D='/home/claude/work/data/'
pt=Pretext(D+'CAE71_hap_merged_assembly_hic_map_mapq0.pretext')
ss=SaveState(D+'CAE71_hap_merged_assembly_hic_map_mapq0.pretext.savestate_1')
N=pt.npix
user_edits=list(ss.edits)
new=[(504,736,+867,False),          # cut H1_scaffold_1 at gap 2, move its 3rd contig behind H2_scaffold_6
     (1837,2188,-925,False),        # H1_scaffold_22 in front of H2_scaffold_6
     (912+317,912+127,0,True),      # invert H1_scaffold_22 px 127..317 in place
     (32130,32172,1264-32130,False)]# H1_scaffold_53 behind H1_scaffold_22
ss.edits=user_edits+new
assert len(ss.edits)%8!=0
ms=MapState(pt)
for e in ss.edits: ms.apply(e)
exp=[('H1_scaffold_1',0,504,0,False),('H2_scaffold_11',504,408,407,True),('H1_scaffold_22',912,127,0,False),('H1_scaffold_22',1039,191,317,True),
     ('H1_scaffold_22',1230,34,318,False),('H1_scaffold_53',1264,43,0,False),('H2_scaffold_6',1307,459,0,False),('H1_scaffold_1',1766,233,504,False),('H2_scaffold_36',1999,233,232,True)]
for i,(nm,s,l,sc,inv) in enumerate(exp):
    c=ms.contigs[i]; assert (pt.names[c[2]],c[0],c[1],c[3],c[4])==(nm,s,l,sc,inv),(i,c)
# scaffold paint: one id per chromosome-scale unit; H1_scaffold_53 left unpainted
ss.scaffs=[(0,1),(1,2),(2,3),(3,3),(4,3),(6,4),(7,5),(8,6)]
HAP1,HAP2,HAPLOTIG=1<<6,1<<7,1<<0
assert ss.tags[6]=='HAP1' and ss.tags[7]=='HAP2' and ss.tags[0]=='Haplotig'
ss.meta_flags=[(0,HAP1),(1,HAP2),(2,HAP1),(3,HAP1),(4,HAP1),(5,HAPLOTIG),(6,HAP2),(7,HAP1),(8,HAP2)]
# waypoints (label shown in PretextView = id+1), placed on the diagonal
wps=[(388,'J'),(834,'J'),(1192,'J'),(1588,'J'),(1039,'inv'),(1230,'inv'),(503,'Y'),(1679,'dup'),(351,'gap1'),(136,'indel')]
wp=b''
for i,(p,_) in enumerate(wps):
    c=p/N-0.5
    wp+=struct.pack('<fffB',c,c,30.0,i)
ss.n_wp=len(wps); ss.waypoints=wp
ss.camera=struct.pack('<fff',1116/N-0.5,0.5-1116/N,13.0)
out='/mnt/user-data/outputs/CAE71_hap_merged_assembly_hic_map_mapq0.pretext.savestate_claude_scaffold1'
payload=ss.dump(out)
print('written',out,os.path.getsize(out),'bytes; payload',len(payload))
np.save(D+'final_pix.npy',ms.pix)
