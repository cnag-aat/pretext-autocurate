"""Align the two haplotypes of the Troglobisium fixture (H2 against H1)."""
import os, sys, time
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate.homology import align_haplotypes, fasta_lengths
from pretext_autocurate import Pretext
fa = os.path.join(R, "seed/fixtures/qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid.fa")
pt = Pretext(os.path.join(R, "seed/fixtures/qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid_mq0.extensions.HR.pretext"), load_graphs=False)
missing = set(pt.names) - set(fasta_lengths(fa))
print("map scaffolds missing from the FASTA:", len(missing), flush=True)
if missing:
    sys.exit(1)
t = time.time()
work = os.path.join(R, "work/qqTroRaco")
paf = align_haplotypes(fa, work, threads=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
for h in ("H1.fa", "H2.fa", "split.done"):          # the split copies are large; the PAF is what is kept
    if os.path.exists(os.path.join(work, h)):
        os.remove(os.path.join(work, h))
print("done", paf, os.path.getsize(paf), "bytes", round(time.time() - t), "s")
