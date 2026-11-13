import sys; sys.path.insert(0,'/home/claude/work/scripts')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def submat(M, idx_r, idx_c=None):
    if idx_c is None: idx_c = idx_r
    return np.asarray(M[np.ix_(idx_r, idx_c)])

def render(sub, out, bounds=None, labels=None, vmax=None, gamma=0.5, figsize=(14,14), dpi=110, extra=None, title=None, vlines=None):
    fig, ax = plt.subplots(figsize=figsize)
    a = sub.astype(np.float32)/255.0
    ax.imshow(a**gamma, cmap='magma_r', vmin=0, vmax=(vmax if vmax else 1.0), interpolation='nearest', origin='upper')
    if bounds is not None:
        for b in bounds:
            ax.axhline(b-0.5, color='k', lw=0.4); ax.axvline(b-0.5, color='k', lw=0.4)
    if vlines:
        for x,c in vlines:
            ax.axvline(x, color=c, lw=0.6, alpha=0.8); ax.axhline(x, color=c, lw=0.6, alpha=0.8)
    if labels:
        ax.set_xticks([p for p,_ in labels]); ax.set_xticklabels([l for _,l in labels], rotation=60, ha='left', fontsize=8)
        ax.xaxis.tick_top()
        ax.set_yticks([p for p,_ in labels]); ax.set_yticklabels([l for _,l in labels], fontsize=8)
    if extra: extra(ax)
    if title: ax.set_title(title, y=-0.04)
    fig.savefig(out, dpi=dpi, bbox_inches='tight'); plt.close(fig)
