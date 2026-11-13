"""Run the haplotype alignment for the bird fixture (long: tens of minutes)."""
import os, sys, time
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate.homology import align_haplotypes
t = time.time()
paf = align_haplotypes(os.path.join(R, "seed/fixtures/bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid.fa"),
                       os.path.join(R, "work/bAegFun"), threads=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
print("done", paf, os.path.getsize(paf), "bytes", round(time.time() - t), "s")
