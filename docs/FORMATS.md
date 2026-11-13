# Binary formats (from sanger-tol/PretextView source; commit in `pretextview_commit.txt`)

All integers little-endian. "deflate" = raw deflate (libdeflate; Python `zlib` with wbits=-15).

## .pretext (PretextMap)

    "pstm" | u32 nComp | u32 nRaw | deflate(header)
    header: u64 total_bp | u32 nContigs | nContigs x { f32 frac ; char name[64] } | u8 texResLog2 | u8 nTexLog2 | u8 nMips
    then nTex*(nTex+1)/2 textures: u32 nBytes | deflate(BC4 data, all mip levels concatenated, level 0 first)
    then extensions, each: "psgh" | u32 nComp | deflate( char name[64] ; i32 value[nPixels] )

- Pixels per side = 2^texResLog2 * 2^nTexLog2 (32768 for the CAE71 maps; 2351 bp/pixel).
- Textures are the upper triangle, row-major: index(x,y) = x*n - x(x-1)/2 + (y-x), x<=y.
  Texture (x,y) holds rows x*res.., columns y*res..; mirror for the lower triangle.
- BC4 (RGTC1): 8 bytes per 4x4 block, 2 endpoints + 16 3-bit indices. Values are 8-bit,
  already normalised by PretextMap: relative intensities, not counts.
- Pixel -> contig: accumulate `frac` **in f32**, pixel_end = (u32)(nPixels * (f64)total). Must be
  replicated exactly or contig boundaries drift by a pixel. Contigs shorter than a pixel vanish.
- Graph names seen: `gap`, `telomere`, `coverage` (also `repeat_density` in other pipelines).

## Savestate

    "ptsx" 0x02 | u8 tail (3 = saved by user, 2 = auto cache) | [u64 headerHash if tail==3] | u32 nComp | u32 nRaw | deflate(payload)

payload, in order:

| bytes | content |
|---|---|
| 2 | theme, setting bits |
| 288 + 28 + 1 + 12 | colours, sizes, colour map id, gamma control points |
| 12 | camera x, y, zoom (f32). y is negated relative to waypoint y |
| 4 | nEdits |
| 12 x nEdits | u32 x, u32 y, i32 delta |
| (nEdits+7)>>3 | invert flags; edit i uses bit (i+1) |
| 1 + 13 x nWayp | nWayp (u8!), then f32 x, f32 y, f32 z, u8 id |
| 4 + 8 x n | painted scaffolds: u32 contigIndex, u32 scaffId |
| 4 + 36 | active tag, tag display data |
| 4 + 12 x n | meta flags: u32 contigIndex, u64 bitmask (bit = tag id) |
| 1 + ... | nTags, then u8 id + NUL-terminated name |
| 32 x nGraphs | per-graph display settings |

Edit semantics (replayed in order on load):
- delta == INT32_MIN: break after pixel x.
- otherwise: move pixel range [x,y] by delta, then if the invert flag is set, reverse the moved range.
  An in-place inversion is delta 0 with the flag set.
- Contig indices for scaffolds/tags are indices **after** all edits, counted along the map. A
  contig is a maximal run of pixels with the same original contig and consecutive coordinates.

Known quirks to handle in the writer:
- When nEdits is a multiple of 8, the last edit's invert bit lives one byte past the flag array
  (it aliases the waypoint-count byte). Avoid by never ending on an inverted edit there, or pad
  with a harmless edit; add a test.
- Waypoint count and ids are u8 (max 255). Edit stack is capped at 32768.
- Waypoint coords: x = y = pixel/nPixels - 0.5 puts it on the diagonal.

## Header hash (ties a savestate to one .pretext file)

    hash = FastHash64(compressedHeaderBytes, seed = FastHash64(filename padded to 64 bytes with NUL, 0xbafd06832de619c2))

It depends on the **file name** and the compressed header, so the MQ0 and MQ30 maps need separate
savestates even though the edits are identical. `seed/scripts/header_hash.py` reproduces the
hash in the curator's own savestate (3ed8d2ab9a753cec for the CAE71 MQ0 map). Renaming a
.pretext file invalidates its savestates.

PretextView can also load an AGP (`Load_AGP`), which is a hash-free way to hand over a result.
