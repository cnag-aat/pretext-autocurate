import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pairs=[('A','B'),('A','C'),('A','D'),('A','E'),('B','C'),('B','D'),('B','E'),('C','D'),('C','E'),('D','E'),('A','A'),('C','C'),('B','B'),('E','E'),('D','D'),('S','A'),('S','B'),('S','C'),('S','D'),('S','E')]
for p,q in pairs:
    X=block(p,q)
    h,w=X.shape
    sc=2 if max(h,w)<800 else 2
    fig,ax=plt.subplots(figsize=(w*sc/100+1.2,h*sc/100+1.0),dpi=100)
    ax.imshow((X/255.0)**0.45,cmap='magma_r',vmin=0,vmax=1,interpolation='nearest',extent=(0,w,h,0),aspect='equal')
    ax.set_xticks(range(0,w,50)); ax.set_yticks(range(0,h,50)); ax.tick_params(labelsize=7)
    ax.grid(True,color='c',lw=0.25,alpha=0.5)
    ax.set_xlabel(f'{q} = {P[q]} (px, native orientation)'); ax.set_ylabel(f'{p} = {P[p]} (px)')
    ax.xaxis.tick_top(); ax.xaxis.set_label_position('top')
    fig.savefig(f'/home/claude/work/out/pair_{p}{q}.png',dpi=100,bbox_inches='tight'); plt.close(fig)
print('ok')
