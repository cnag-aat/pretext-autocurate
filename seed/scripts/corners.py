import sys; sys.path.insert(0,'/home/claude/work/scripts')
from junction2 import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def corner_grid(name, arms, out, smin=1, smax=60):
    keys=list(arms.keys()); n=len(keys)
    fig,axs=plt.subplots(n,n,figsize=(3.2*n,3.2*n))
    for i,x in enumerate(keys):
        for j,y in enumerate(keys):
            ax=axs[i,j]
            if i==j: ax.axis('off'); continue
            T=arm_contact(arms[x],arms[y],smin,smax)
            ax.imshow((np.clip(T,0,400)/400)**0.5,cmap='magma_r',vmin=0,vmax=1,origin='lower',extent=(smin,smax,smin,smax))
            ax.set_title(f'{x} (y) x {y} (x)',fontsize=9); ax.tick_params(labelsize=6)
    fig.suptitle(name+': collapsed contacts vs distance from junction (px)')
    fig.savefig(out,dpi=90,bbox_inches='tight'); plt.close(fig)
corner_grid('Junction J',J,'/home/claude/work/out/cornerJ.png')
corner_grid('Junction K',K,'/home/claude/work/out/cornerK.png')
