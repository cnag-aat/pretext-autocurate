"""miniprot of the aves_odb12 BUSCO proteins against the bird assembly (long)."""
import os, sys, time
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate import busco
t = time.time()
d = busco.resolve_lineage(sys.argv[2] if len(sys.argv) > 2 else "aves_odb12")
out = busco.run_miniprot(os.path.join(R, "seed/fixtures/bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid.fa"),
                         busco.lineage_proteins(d), os.path.join(R, "work/bAegFun/aves_odb12.ancestral.miniprot.paf"),
                         threads=int(sys.argv[1]) if len(sys.argv) > 1 else 6)
print("done", out, os.path.getsize(out), "bytes", round(time.time() - t), "s")
