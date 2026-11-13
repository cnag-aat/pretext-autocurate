"""Savestate header hash, as PretextView computes it in LoadFile (verified on CAE71)."""
import struct, sys, os
M64 = (1 << 64) - 1
m = 0x880355f21e6d1965

def _mix(h):
    h ^= h >> 23; h = (h * 0x2127599bf4325c37) & M64; h ^= h >> 47
    return h

def fasthash64(buf, seed):
    n = len(buf); h = seed ^ ((n * m) & M64)
    for i in range(n // 8):
        v, = struct.unpack_from('<Q', buf, 8 * i)
        h ^= _mix(v); h = (h * m) & M64
    tail = buf[8 * (n // 8):]
    if tail:
        v = int.from_bytes(tail, 'little')
        h ^= _mix(v); h = (h * m) & M64
    return _mix(h)

def header_hash(pretext_path):
    name = os.path.basename(pretext_path).encode()[:64].ljust(64, b'\0')
    with open(pretext_path, 'rb') as f:
        assert f.read(4) == b'pstm'
        nc, _ = struct.unpack('<II', f.read(8)); comp = f.read(nc)
    return fasthash64(comp, fasthash64(name, 0xbafd06832de619c2))

if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(struct.pack('<Q', header_hash(p)).hex(), p)
