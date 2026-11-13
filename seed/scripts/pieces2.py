import sys; sys.path.insert(0,'/home/claude/work/scripts')
from pieces import *
import homology
P.update({'F':'H1_scaffold_23','G':'H2_scaffold_28','H':'H2_scaffold_29','I':'H2_scaffold_51'})
for k in 'FGHI': R[k]=crange(P[k])
