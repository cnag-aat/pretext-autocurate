import sys, time; sys.path.insert(0,'/home/claude/work/scripts')
from pretext import *
import numpy as np
src, out = sys.argv[1], sys.argv[2]
pt=Pretext(src)
N=pt.npix; R=pt.tex_res
M=np.lib.format.open_memmap(out, mode='w+', dtype=np.uint8, shape=(N,N))
t=time.time()
for x in range(pt.n_tex):
    for y in range(x, pt.n_tex):
        tex=pt.read_tex(x,y)
        M[x*R:(x+1)*R, y*R:(y+1)*R]=tex
        if x!=y: M[y*R:(y+1)*R, x*R:(x+1)*R]=tex.T
M.flush()
print('done',time.time()-t, pt.n_contigs, pt.total_len)
