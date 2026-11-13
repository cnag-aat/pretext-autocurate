# Test case 1: CAE71 hap-merged, "scaffold 1" group (hand-derived, session 1)

Map: `CAE71_hap_merged_assembly_hic_map_mapq0.pretext` (147 scaffolds, 77.05 Mb, hifiasm --ont
hap1+hap2 concatenated, mixed non-clonal protist culture). Coordinates below are pixels in each
scaffold's native orientation (2351 bp/pixel).

Pieces: A=H1_scaffold_1 (737 px, gaps at 351 and 503), B=H2_scaffold_11 (408), C=H2_scaffold_6 (459),
D=H2_scaffold_36 (233), E=H1_scaffold_22 (352), S=H1_scaffold_53 (43).

Homology (MQ0 off-diagonal traces):
- A[0-388] ~ B[383-76] reversed
- A[386-479] ~ C[280-371]
- C[22-135] ~ E[0-125]; C[137-281] ~ E[316-163] reversed
- B[16-78] ~ E[124-165]; B[0-14] ~ E[330-318] reversed; B[0-22] ~ S[22-0]
- A[504-736] ~ D[232-1] reversed
- A[502-545] ~ C[458-418] reversed; A[544-568] ~ C[382-402]; repeat A[482-490] ~ A[561-569]

Findings the tool should reproduce:
1. Four-way junction at A388 / B77 / C280 / E165: arms P,Q,U,V with joins P-Q (A), P-U (B),
   V-Q (C), V-U (E) all at cis level (121-187) and P-V, Q-U at background (8-19).
2. E[127-317] inverted relative to Hi-C: as assembled its breakpoints score 53 and 17, flipped
   they score 169 and 111. Edit: invert in place.
3. A is misjoined at gap 2 (503|504): contact between A[..480] and A[569..] is background (7.8).
   A[504-736] is the homolog of D. Edit: cut and move next to D.
4. A gap 1 (351) is fine (B is contiguous across the homologous position).
5. A ~135 kb block (A[504-561]) is present three times: start of A[504-736], start of D, and
   rearranged at the telomere end of C.
6. The end of A before gap 2 has cis-level contact with H2_scaffold_29 and H1_scaffold_23
   (not resolved in session 1).
7. S is a low-confidence placement after E.

Expected first contigs after the edits (start px, length, scaffold, start coord, inverted):

    0    504  H1_scaffold_1   0    +
    504  408  H2_scaffold_11  407  -
    912  127  H1_scaffold_22  0    +
    1039 191  H1_scaffold_22  317  -
    1230 34   H1_scaffold_22  318  +
    1264 43   H1_scaffold_53  0    +
    1307 459  H2_scaffold_6   0    +
    1766 233  H1_scaffold_1   504  +
    1999 233  H2_scaffold_36  232  -

Fixtures: `seed/fixtures/*.savestate_1` (curator's 97 edits), `*.savestate_claude_scaffold1`
(those plus 4 edits, paint, tags, 10 waypoints) for MQ0 and MQ30. The MQ0 one has not yet been
opened in PretextView by anyone; confirming that it loads is the first task.

Caveat: these findings came from one analysis pass and have not been checked against sequence
alignments. Treat them as a regression target for the method, not as biological truth.
