import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pretext import *
import numpy as np
PT=Pretext('/home/claude/work/data/CAE71_hap_merged_assembly_hic_map_mapq0.pretext')
M=np.load('/home/claude/work/data/mq0_full.npy', mmap_mode='r')
M30=np.load('/home/claude/work/data/mq30_full.npy', mmap_mode='r')
def crange(name):
    c=PT.names.index(name); w=np.where(PT.orig_contig==c)[0]; return int(w[0]), int(w[-1])+1
P={'A':'H1_scaffold_1','B':'H2_scaffold_11','C':'H2_scaffold_6','D':'H2_scaffold_36','E':'H1_scaffold_22','S':'H1_scaffold_53'}
R={k:crange(v) for k,v in P.items()}
def block(p,q,mat=None):
    mat = M if mat is None else mat
    (a,b),(c,d)=R[p],R[q]
    return np.asarray(mat[a:b,c:d]).astype(np.float32)
