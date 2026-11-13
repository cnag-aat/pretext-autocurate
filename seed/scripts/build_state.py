import sys, struct; sys.path.insert(0,'/home/claude/work/scripts')
from pretext import *
import numpy as np
D='/home/claude/work/data/'
pt=Pretext(D+'CAE71_hap_merged_assembly_hic_map_mapq0.pretext')
ss=SaveState(D+'CAE71_hap_merged_assembly_hic_map_mapq0.pretext.savestate_1')
# --- round-trip check of the writer on the untouched state
raw0=ss.raw
import io, os
tmp='/tmp/claude-0/-home-claude/9718f7b1-9961-575b-a492-e45b9211bdf2/scratchpad/rt.savestate'
os.makedirs(os.path.dirname(tmp),exist_ok=True)
raw1=ss.dump(tmp)
assert raw1==raw0, 'round trip mismatch'
print('writer round-trip on your savestate: identical payload (%d bytes)'%len(raw1))

ms=MapState(pt)
for e in ss.edits: ms.apply(e)
def where(name):
    c=pt.names.index(name)
    return [(i,s,l,inv,sc) for i,(s,l,ob,sc,inv,o) in enumerate(ms.contigs) if ob==c]
for n in ['H1_scaffold_1','H2_scaffold_11','H2_scaffold_6','H2_scaffold_36','H1_scaffold_22','H1_scaffold_53']:
    print(n, where(n))
new=[]
# 1. cut H1_scaffold_1 at gap 2 (px 503|504) and move the 3rd contig (504-736) behind H2_scaffold_6
new.append((504,736,+867,False))
# 2. move H1_scaffold_22 (now 1837-2188) in front of H2_scaffold_6
new.append((1837,2188,-925,False))
# 3. invert H1_scaffold_22 local px 127..317 in place  (stored as max,min like PretextView does)
new.append((912+317,912+127,0,True))
# 4. bring H1_scaffold_53 (32130-32172) in behind H1_scaffold_22
new.append((32130,32172,1264-32130,False))
for e in new: ms.apply(e)
print('--- first contigs after my edits')
for i,(s,l,ob,sc,inv,o) in enumerate(ms.contigs[:11]):
    print(i,s,s+l-1,l,pt.names[ob],'rel',sc,'-' if inv else '+')
np.save(D+'final_pix.npy',ms.pix)
import pickle; pickle.dump((new,ms.contigs[:12]),open(D+'final_edits.pkl','wb'))
