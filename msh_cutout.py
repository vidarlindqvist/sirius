# msh_cutout.py
#
# Cuts a box-shaped region out of a large ASCII Fluent .msh file and writes it
# as VTK files that ParaView opens directly. It parses the .msh itself, so
# ParaView's Fluent reader is never involved. Needs only numpy.
#
# Run with ParaView's interpreter (numpy is bundled):
#   & "C:\Users\vidlin-2\AppData\Local\Programs\ParaView\ParaView-6.2.0-Windows-Python3.12-msvc2017-AMD64\bin\pvpython.exe" C:\temp\msh_cutout.py
#
# Optional command line, overrides the settings below:
#   msh_cutout.py [--info] [--inside] MESH OUT_BASE  XMIN YMIN ZMIN  XMAX YMAX ZMAX
#
# Output: OUT_BASE.vtm  (open this one)  +  folder OUT_BASE\ with one .vtu per zone
#   volume/<cell zone>     tets, prisms, pyramids, hexes
#   surfaces/<face zone>   boundary faces (walls, symmetry, ...) of the kept cells
#   Cell array NodesPerElement: 3 tri, 4 quad | 4 tet, 5 pyramid, 6 prism, 8 hex

import os
import re
import sys
import time

import numpy as np

# ---------------------------------------------------------------- settings
MESH_FILE = r'C:\temp\fHSP_QUALITY_20.msh'
OUT_BASE  = r'C:\temp\wing_region'

# Box as two opposite corners, in mesh units. Order does not matter.
BOX_CORNER_A = [-1.0, -2.0, -5.0]
BOX_CORNER_B = [19.0,  3.0,  0.0]

# 'touch'  : keep every cell with at least one node inside the box
# 'inside' : keep only cells lying completely inside the box
MODE = 'touch'

# True: read the file and print the zone table only, write nothing.
INFO_ONLY = False
# -------------------------------------------------------------------------

BLOCK = 32 * 1024 * 1024          # bytes read per step

VTK_TRI, VTK_QUAD, VTK_TET, VTK_HEX, VTK_WEDGE, VTK_PYRAMID = 5, 9, 10, 12, 13, 14
BC_INTERIOR, BC_PARENT = 2, 31
BC_NAMES = {2: 'interior', 3: 'wall', 4: 'pressure-inlet', 5: 'pressure-outlet',
            7: 'symmetry', 8: 'periodic-shadow', 9: 'pressure-far-field',
            10: 'velocity-inlet', 12: 'periodic', 14: 'fan', 20: 'mass-flow-inlet',
            24: 'interface', 31: 'parent', 36: 'outflow', 37: 'axis'}

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.0f} s] {msg}", flush=True)


# ============================================================ file scanning
_NONWS = re.compile(rb'\S')
_PQ = re.compile(rb'[()"]')
_SECID = re.compile(rb'\(\s*(\d+)')
_ZONE = re.compile(rb'\(\s*(?:39|45)\s*\(\s*(\d+)\s+([^\s()]+)\s+([^\s()]+)')


class Stream:
    """Forward-only reader over a huge text file, a block at a time."""

    def __init__(self, path):
        self.f = open(path, 'rb', buffering=0)
        self.size = os.path.getsize(path)
        self.buf = b''
        self.pos = 0            # read position in buf
        self.base = 0           # file offset of buf[0]
        self.eof = False

    def tell(self):
        return self.base + self.pos

    def _fill(self):
        if self.eof:
            return False
        data = self.f.read(BLOCK)
        if not data:
            self.eof = True
            return False
        self.base += self.pos
        self.buf = self.buf[self.pos:] + data
        self.pos = 0
        return True

    def need(self, n):
        while len(self.buf) - self.pos < n and self._fill():
            pass

    def skip_ws(self):
        """Skip whitespace; return the next byte (b'' at end of file)."""
        while True:
            m = _NONWS.search(self.buf, self.pos)
            if m:
                self.pos = m.start()
                return self.buf[self.pos:self.pos + 1]
            self.pos = len(self.buf)
            if not self._fill():
                return b''

    def read_group(self, keep=True):
        """Consume one balanced (...) group starting at pos; return its text."""
        depth, in_str, scan = 0, False, self.pos
        while True:
            for m in _PQ.finditer(self.buf, scan):
                c = m.group()
                if c == b'"':
                    in_str = not in_str
                elif in_str:
                    continue
                elif c == b'(':
                    depth += 1
                else:
                    depth -= 1
                    if depth == 0:
                        text = self.buf[self.pos:m.end()] if keep else b''
                        self.pos = m.end()
                        return text
            if keep:
                scan = len(self.buf) - self.pos
            else:
                self.pos = len(self.buf)
                scan = 0
            if not self._fill():
                raise ValueError("file ends in the middle of a section (truncated?)")

    def body_chunks(self):
        """Yield the data up to the next ')' in pieces that end on a line break."""
        while True:
            end = self.buf.find(b')', self.pos)
            if end >= 0:
                if end > self.pos:
                    yield self.buf[self.pos:end]
                self.pos = end + 1
                return
            cut = self.buf.rfind(b'\n', self.pos)
            if cut < 0:
                cut = self.buf.rfind(b' ', self.pos)
            if cut >= self.pos:
                yield self.buf[self.pos:cut + 1]
                self.pos = cut + 1
            if not self._fill():
                raise ValueError("file ends in the middle of a data section (truncated?)")


# ============================================================ number parsing
_HEXVAL = np.full(256, -2, np.int8)             # -2 invalid, -1 whitespace
for _c in b' \t\r\n':
    _HEXVAL[_c] = -1
for _i, _c in enumerate(b'0123456789'):
    _HEXVAL[_c] = _i
for _i, _c in enumerate(b'abcdef'):
    _HEXVAL[_c] = 10 + _i
for _i, _c in enumerate(b'ABCDEF'):
    _HEXVAL[_c] = 10 + _i


def parse_hex(chunk):
    """All hexadecimal integers in a bytes chunk -> (values int64, start offsets, raw bytes)."""
    a = np.frombuffer(chunk, dtype=np.uint8)
    d = _HEXVAL[a]
    if (d == -2).any():
        bad = int(np.flatnonzero(d == -2)[0])
        raise ValueError(f"unexpected character in integer data: {bytes(chunk[max(0, bad - 30):bad + 30])!r}")
    isd = (d >= 0).view(np.int8)
    edge = np.diff(isd, prepend=np.int8(0), append=np.int8(0))
    starts = np.flatnonzero(edge == 1)
    lens = np.flatnonzero(edge == -1) - starts
    val = np.zeros(len(starts), np.int64)
    for k in range(int(lens.max()) if len(lens) else 0):
        m = lens > k
        if m.all():
            val = val * 16 + d[starts + k]
        else:
            val[m] = val[m] * 16 + d[starts[m] + k]
    return val, starts, a


def parse_floats(chunk):
    try:
        v = np.fromstring(chunk.decode('ascii'), dtype=np.float64, sep=' ')
    except Exception:
        v = None
    if v is None or len(v) == 0:
        v = np.array(chunk.split(), dtype=np.float64)
    return v


# ============================================================ reading the mesh
def read_msh(path):
    s = Stream(path)
    xyz = None
    n_cells = 0
    cell_zones = []          # (zone, first, last)
    face_blocks = []         # dict(zone, bc, k, rows)   rows = [node0..node(k-1), c0, c1]
    zone_names = {}          # id -> (type, name)
    bad_cells = []           # cells touching faces this script cannot use (polygons)
    n_poly = 0
    next_report = 0.1

    def progress():
        nonlocal next_report
        frac = s.tell() / max(s.size, 1)
        if frac >= next_report:
            log(f"  {100 * frac:3.0f} % of file read")
            next_report = (int(frac * 10) + 1) / 10

    while True:
        c = s.skip_ws()
        if not c:
            break
        if c != b'(':
            raise ValueError(f"unexpected text at byte {s.tell()}: {s.buf[s.pos:s.pos + 40]!r}")
        s.need(64)
        m = _SECID.match(s.buf, s.pos)
        sec = int(m.group(1)) if m else -1

        if sec in (10, 12, 13):
            s.pos = m.end()
            s.skip_ws()
            hdr = [int(t, 16) for t in s.read_group()[1:-1].split()]
            has_body = s.skip_ws() == b'('
            if has_body:
                s.pos += 1
            zone, first, last = hdr[0], hdr[1], hdr[2]

            if sec == 10:
                if len(hdr) > 4 and hdr[4] != 3:
                    raise ValueError("this is a 2D mesh; only 3D meshes are supported")
                if xyz is None or last > len(xyz):
                    new = np.full((last, 3), np.nan)
                    if xyz is not None:
                        new[:len(xyz)] = xyz
                    xyz = new
                if has_body:
                    log(f"reading {last - first + 1:,} nodes ...")
                    flat, at = xyz.reshape(-1), (first - 1) * 3
                    for chunk in s.body_chunks():
                        v = parse_floats(chunk)
                        flat[at:at + len(v)] = v
                        at += len(v)
                        progress()
                    if at != last * 3:
                        raise ValueError(f"node section: expected {(last - first + 1) * 3} numbers, "
                                         f"got {at - (first - 1) * 3}")

            elif sec == 12:
                if zone == 0:
                    n_cells = max(n_cells, last)
                else:
                    cell_zones.append((zone, first, last))
                    n_cells = max(n_cells, last)
                if has_body:                       # per-cell element types, not needed
                    for _ in s.body_chunks():
                        pass

            else:  # 13, faces
                bc = hdr[3] if len(hdr) > 3 else 0
                ft = hdr[4] if len(hdr) > 4 else 0
                if has_body:
                    if zone and not face_blocks:
                        log("reading faces ...")
                    rows = {3: [], 4: []}
                    if ft in (3, 4):
                        w, rem = ft + 2, np.zeros(0, np.int64)
                        for chunk in s.body_chunks():
                            t = parse_hex(chunk)[0]
                            if len(rem):
                                t = np.concatenate((rem, t))
                            n = len(t) // w * w
                            rem = t[n:]
                            rows[ft].append(t[:n].reshape(-1, w).astype(np.int32))
                            progress()
                        if len(rem):
                            raise ValueError(f"face zone {zone}: incomplete last face")
                    elif ft in (0, 5):             # count-prefixed, one face per line
                        for chunk in s.body_chunks():
                            t, starts, a = parse_hex(chunk)
                            if not len(t):
                                continue
                            line = np.searchsorted(np.flatnonzero(a == 10), starts)
                            f0 = np.flatnonzero(np.diff(line, prepend=line[0] - 1))
                            cnt = t[f0]
                            if not (np.diff(f0, append=len(t)) == cnt + 3).all():
                                raise ValueError(f"face zone {zone}: unexpected line layout "
                                                 f"in mixed-face section")
                            for k in (3, 4):
                                f = f0[cnt == k]
                                if len(f):
                                    rows[k].append(t[f[:, None] + np.arange(1, k + 3)].astype(np.int32))
                            f = f0[(cnt != 3) & (cnt != 4)]
                            if len(f):
                                n_poly += len(f)
                                bad_cells.append(t[f + cnt[(cnt != 3) & (cnt != 4)] + 1])
                                bad_cells.append(t[f + cnt[(cnt != 3) & (cnt != 4)] + 2])
                            progress()
                    else:
                        raise ValueError(f"face zone {zone}: face type {ft} not supported "
                                         f"(2D mesh?)")
                    for k in (3, 4):
                        if rows[k]:
                            r = np.concatenate(rows[k]) if len(rows[k]) > 1 else rows[k][0]
                            r[:, :k] -= 1                       # nodes 0-based
                            face_blocks.append(dict(zone=zone, bc=bc, k=k, rows=r))
                    del rows

            if has_body:
                if s.skip_ws() != b')':
                    raise ValueError(f"section {sec}: missing closing parenthesis at byte {s.tell()}")
                s.pos += 1
            else:
                if s.skip_ws() != b')':
                    raise ValueError(f"section {sec}: malformed header at byte {s.tell()}")
                s.pos += 1

        elif sec in (39, 45):
            text = s.read_group()
            zm = _ZONE.match(text)
            if zm:
                zone_names[int(zm.group(1))] = (zm.group(2).decode('ascii', 'replace'),
                                                zm.group(3).decode('ascii', 'replace'))
        elif sec in (2010, 3010, 2012, 3012, 2013, 3013, 2041, 3041):
            raise ValueError("this .msh is written in binary; this script reads ASCII only. "
                             "Re-export the mesh as ASCII.")
        else:
            s.read_group(keep=False)
        progress()

    s.f.close()
    if xyz is None or not face_blocks:
        raise ValueError("no nodes or no faces found; is this a Fluent mesh file?")
    if np.isnan(xyz).any():
        raise ValueError("some node coordinates are missing from the file")
    for b in face_blocks:
        n_cells = max(n_cells, int(b['rows'][:, b['k']:].max()))
    bad = np.unique(np.concatenate(bad_cells)) if bad_cells else np.zeros(0, np.int64)
    return xyz, n_cells, cell_zones, face_blocks, zone_names, bad, n_poly


# ============================================================ cell rebuilding
def det3(p, a, b, c, d):
    """Signed 6*volume of the tetrahedron spanned by b-a, c-a, d-a."""
    pa = p[a]
    return np.einsum('ij,ij->i', np.cross(p[b] - pa, p[c] - pa), p[d] - pa)


def opposite_nodes(A, Q):
    """A (n,m): nodes of one face. Q (n,q,4): the cell's other quad faces.
    For every node of A, the node it is joined to by an edge leaving A."""
    in_a = (Q[:, :, :, None] == A[:, None, None, :]).any(axis=3)
    cand = np.where(np.roll(in_a, -1, axis=2), np.roll(Q, 1, axis=2), np.roll(Q, -1, axis=2))
    B = np.empty_like(A)
    ok = np.ones(len(A), bool)
    for i in range(A.shape[1]):
        hit = Q == A[:, i][:, None, None]
        good = hit.reshape(len(A), -1).sum(axis=1) == 2
        ok &= good
        hit[~good] = False
        hit[~good, 0, :2] = True                   # placeholder so the reshape works
        B[:, i] = cand[hit].reshape(-1, 2)[:, 0]
    return B, ok


def build_cells(xyz, tri_c, tri_n, quad_c, quad_n, n_kept, bad_kept):
    """Rebuild VTK cells from their faces.
    tri_c/quad_c: compact cell id per face entry, tri_n/quad_n: face nodes.
    Returns list of (vtk_type, conn, cell ids) and the number of unusable cells."""
    o = np.argsort(tri_c, kind='stable')
    tri_c, tri_n = tri_c[o], tri_n[o]
    o = np.argsort(quad_c, kind='stable')
    quad_c, quad_n = quad_c[o], quad_n[o]
    nt = np.bincount(tri_c, minlength=n_kept)
    nq = np.bincount(quad_c, minlength=n_kept)
    usable = ~bad_kept
    out, n_done = [], 0

    def faces_of(cls, c, n, per):
        sel = cls[c]
        return n[sel].reshape(-1, per, n.shape[1]), c[sel][::per]

    # tetrahedra: 4 triangles
    cls = usable & (nt == 4) & (nq == 0)
    if cls.any():
        T, ids = faces_of(cls, tri_c, tri_n, 4)
        base, other = T[:, 0], T[:, 1]
        free = (other[:, :, None] != base[:, None, :]).all(axis=2)
        ok = free.sum(axis=1) == 1
        T, ids, base, other, free = T[ok], ids[ok], base[ok], other[ok], free[ok]
        conn = np.column_stack((base, other[free]))
        flip = det3(xyz, conn[:, 0], conn[:, 1], conn[:, 2], conn[:, 3]) < 0
        conn[flip] = conn[flip][:, [0, 2, 1, 3]]
        out.append((VTK_TET, conn, ids))
        n_done += len(ids)

    # pyramids: 1 quad + 4 triangles
    cls = usable & (nt == 4) & (nq == 1)
    if cls.any():
        T, ids = faces_of(cls, tri_c, tri_n, 4)
        Q, ids2 = faces_of(cls, quad_c, quad_n, 1)
        assert (ids == ids2).all()
        base, tri = Q[:, 0], T[:, 0]
        free = (tri[:, :, None] != base[:, None, :]).all(axis=2)
        ok = free.sum(axis=1) == 1
        ids, base, tri, free = ids[ok], base[ok], tri[ok], free[ok]
        conn = np.column_stack((base, tri[free]))
        flip = det3(xyz, conn[:, 0], conn[:, 1], conn[:, 3], conn[:, 4]) < 0
        conn[flip] = conn[flip][:, [0, 3, 2, 1, 4]]
        out.append((VTK_PYRAMID, conn, ids))
        n_done += len(ids)

    # prisms (wedges): 2 triangles + 3 quads
    cls = usable & (nt == 2) & (nq == 3)
    if cls.any():
        T, ids = faces_of(cls, tri_c, tri_n, 2)
        Q, ids2 = faces_of(cls, quad_c, quad_n, 3)
        assert (ids == ids2).all()
        A = T[:, 0].copy()
        B, ok = opposite_nodes(A, Q)
        ids, A, B = ids[ok], A[ok], B[ok]
        conn = np.column_stack((A, B))
        # VTK wants the normal of triangle (0,1,2) pointing towards (3,4,5)
        flip = det3(xyz, conn[:, 0], conn[:, 1], conn[:, 2], conn[:, 3]) < 0
        conn[flip] = conn[flip][:, [0, 2, 1, 3, 5, 4]]
        out.append((VTK_WEDGE, conn, ids))
        n_done += len(ids)

    # hexahedra: 6 quads
    cls = usable & (nt == 0) & (nq == 6)
    if cls.any():
        Q, ids = faces_of(cls, quad_c, quad_n, 6)
        A = Q[:, 0].copy()
        B, ok = opposite_nodes(A, Q[:, 1:])
        ids, A, B = ids[ok], A[ok], B[ok]
        conn = np.column_stack((A, B))
        flip = det3(xyz, conn[:, 0], conn[:, 1], conn[:, 3], conn[:, 4]) < 0
        conn[flip] = conn[flip][:, [0, 3, 2, 1, 4, 7, 6, 5]]
        out.append((VTK_HEX, conn, ids))
        n_done += len(ids)

    return out, n_kept - n_done


# ============================================================ VTK writers
def write_vtu(path, xyz, pieces):
    """pieces: list of (vtk_type, conn (n,k) with 0-based global node ids)."""
    pieces = [(t, c) for t, c in pieces if len(c)]
    flat = np.concatenate([c.reshape(-1) for _, c in pieces])
    used = np.zeros(len(xyz), bool)
    used[flat] = True
    idx = np.flatnonzero(used)
    lut = np.empty(len(xyz), np.int64)
    lut[idx] = np.arange(len(idx))
    per = np.concatenate([np.full(len(c), c.shape[1], np.int64) for _, c in pieces])
    arrays = [
        ('Float64', np.ascontiguousarray(xyz[idx], '<f8')),
        ('Int64', lut[flat].astype('<i8')),
        ('Int64', np.cumsum(per).astype('<i8')),
        ('UInt8', np.concatenate([np.full(len(c), t, np.uint8) for t, c in pieces])),
        ('Int32', per.astype('<i4')),
    ]
    offs, at = [], 0
    for _, arr in arrays:
        offs.append(at)
        at += 8 + arr.nbytes
    xml = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n'
        ' <UnstructuredGrid>\n'
        f'  <Piece NumberOfPoints="{len(idx)}" NumberOfCells="{len(per)}">\n'
        '   <Points>\n'
        f'    <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="appended" offset="{offs[0]}"/>\n'
        '   </Points>\n'
        '   <Cells>\n'
        f'    <DataArray type="Int64" Name="connectivity" format="appended" offset="{offs[1]}"/>\n'
        f'    <DataArray type="Int64" Name="offsets" format="appended" offset="{offs[2]}"/>\n'
        f'    <DataArray type="UInt8" Name="types" format="appended" offset="{offs[3]}"/>\n'
        '   </Cells>\n'
        '   <CellData>\n'
        f'    <DataArray type="Int32" Name="NodesPerElement" format="appended" offset="{offs[4]}"/>\n'
        '   </CellData>\n'
        '  </Piece>\n'
        ' </UnstructuredGrid>\n'
        ' <AppendedData encoding="raw">\n_'
    )
    with open(path, 'wb') as f:
        f.write(xml.encode('ascii'))
        for _, arr in arrays:
            f.write(np.uint64(arr.nbytes).tobytes())
            arr.tofile(f)
        f.write(b'\n </AppendedData>\n</VTKFile>\n')
    return len(idx), len(per)


def write_vtm(path, groups):
    """groups: list of (group name, [(block name, relative file)])."""
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">',
             ' <vtkMultiBlockDataSet>']
    for gi, (gname, items) in enumerate(groups):
        lines.append(f'  <Block index="{gi}" name="{gname}">')
        for i, (name, rel) in enumerate(items):
            lines.append(f'   <DataSet index="{i}" name="{name}" file="{rel}"/>')
        lines.append('  </Block>')
    lines += [' </vtkMultiBlockDataSet>', '</VTKFile>', '']
    with open(path, 'w', encoding='ascii') as f:
        f.write('\n'.join(lines))


def safe(name):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)


# ============================================================ main
def main():
    mesh, out_base, mode, info_only = MESH_FILE, OUT_BASE, MODE, INFO_ONLY
    ca, cb = BOX_CORNER_A, BOX_CORNER_B
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--info' in sys.argv:
        info_only = True
    if '--inside' in sys.argv:
        mode = 'inside'
    if len(args) >= 1:
        mesh = args[0]
    if len(args) >= 2:
        out_base = args[1]
    if len(args) >= 8:
        ca, cb = [float(v) for v in args[2:5]], [float(v) for v in args[5:8]]
    if mode not in ('touch', 'inside'):
        raise SystemExit("MODE must be 'touch' or 'inside'")
    lo = np.minimum(ca, cb).astype(float)
    hi = np.maximum(ca, cb).astype(float)
    if (hi - lo).min() <= 0:
        raise SystemExit("Box has zero thickness in at least one direction.")

    log(f"file: {mesh}  ({os.path.getsize(mesh) / 1e9:.2f} GB)")
    xyz, n_cells, cell_zones, blocks, names, bad, n_poly = read_msh(mesh)
    n_faces = sum(len(b['rows']) for b in blocks)
    log(f"read done: {len(xyz):,} nodes, {n_faces:,} faces, {n_cells:,} cells")
    dmin, dmax = xyz.min(axis=0), xyz.max(axis=0)
    print(f"\n  domain  x {dmin[0]:.4g} .. {dmax[0]:.4g}   y {dmin[1]:.4g} .. {dmax[1]:.4g}"
          f"   z {dmin[2]:.4g} .. {dmax[2]:.4g}")
    print(f"  box     x {lo[0]:.4g} .. {hi[0]:.4g}   y {lo[1]:.4g} .. {hi[1]:.4g}"
          f"   z {lo[2]:.4g} .. {hi[2]:.4g}   (mode: {mode})\n")
    if n_poly:
        print(f"  note: {n_poly:,} polygonal faces found; cells using them are skipped\n")

    def zname(z):
        return names[z][1] if z in names else f"zone_{z}"

    def ztype(z, bc):
        return names[z][0] if z in names else BC_NAMES.get(bc, f"type_{bc}")

    # ---- which cells to keep
    inbox = ((xyz >= lo) & (xyz <= hi)).all(axis=1)
    keep = np.zeros(n_cells + 1, bool)
    if mode == 'touch':
        for b in blocks:
            k, r = b['k'], b['rows']
            hit = inbox[r[:, :k]].any(axis=1)
            keep[r[hit, k]] = True
            keep[r[hit, k + 1]] = True
    else:
        seen = np.zeros(n_cells + 1, bool)
        outside = np.zeros(n_cells + 1, bool)
        for b in blocks:
            k, r = b['k'], b['rows']
            seen[r[:, k]] = True
            seen[r[:, k + 1]] = True
            miss = ~inbox[r[:, :k]].all(axis=1)
            outside[r[miss, k]] = True
            outside[r[miss, k + 1]] = True
        keep = seen & ~outside
    keep[0] = False
    n_kept = int(keep.sum())
    lut = np.full(n_cells + 1, -1, np.int32)
    lut[keep] = np.arange(n_kept, dtype=np.int32)

    # ---- collect faces of kept cells, and boundary faces per zone
    ent = {3: ([], []), 4: ([], [])}
    surf = {}                      # zone -> {k: [nodes]}
    table = {}                     # zone -> [bc, n faces, n faces kept, min, max]
    for b in blocks:
        k, r, z = b['k'], b['rows'], b['zone']
        nodes = r[:, :k]
        boundary = b['bc'] not in (BC_INTERIOR, BC_PARENT)
        m_any = np.zeros(len(r), bool)
        for side in (k, k + 1):
            c = r[:, side]
            m = keep[c]
            m_any |= m
            if m.any():
                ent[k][0].append(lut[c[m]])
                ent[k][1].append(nodes[m])
        row = table.setdefault(z, [b['bc'], 0, 0, None, None])
        row[1] += len(r)
        row[2] += int(m_any.sum())
        if boundary:
            p = xyz[np.unique(nodes)]
            mn, mx = p.min(axis=0), p.max(axis=0)
            row[3] = mn if row[3] is None else np.minimum(row[3], mn)
            row[4] = mx if row[4] is None else np.maximum(row[4], mx)
            if m_any.any():
                surf.setdefault(z, {}).setdefault(k, []).append(nodes[m_any])

    print("  face zones (extents are for the whole zone, not just the box):")
    print(f"  {'id':>5}  {'name':<28} {'type':<18} {'faces':>12} {'in cut-out':>11}   extents")
    for z in sorted(table):
        bc, n, nk, mn, mx = table[z]
        ext = ('' if mn is None else
               f"x {mn[0]:.4g}..{mx[0]:.4g}  y {mn[1]:.4g}..{mx[1]:.4g}  z {mn[2]:.4g}..{mx[2]:.4g}")
        print(f"  {z:>5}  {zname(z):<28.28} {ztype(z, bc):<18.18} {n:>12,} {nk:>11,}   {ext}")
    print()

    if info_only:
        log(f"info only: {n_kept:,} cells would be kept; nothing written")
        return
    if n_kept == 0:
        raise SystemExit("No cells in the box. Compare the box with the domain and zone "
                         "extents printed above and adjust the corners.")

    # ---- rebuild cells
    log(f"rebuilding {n_kept:,} cells ...")
    cat = lambda lst, shape, dt: (np.concatenate(lst) if lst else np.zeros(shape, dt))
    tri_c, tri_n = cat(ent[3][0], 0, np.int32), cat(ent[3][1], (0, 3), np.int32)
    quad_c, quad_n = cat(ent[4][0], 0, np.int32), cat(ent[4][1], (0, 4), np.int32)
    del ent
    bad_kept = np.zeros(n_kept, bool)
    if len(bad):
        b_ = lut[bad[(bad > 0) & (bad <= n_cells)]]
        bad_kept[b_[b_ >= 0]] = True
    cells, n_skipped = build_cells(xyz, tri_c, tri_n, quad_c, quad_n, n_kept, bad_kept)
    del tri_c, tri_n, quad_c, quad_n
    label = {VTK_TET: 'tets', VTK_PYRAMID: 'pyramids', VTK_WEDGE: 'prisms', VTK_HEX: 'hexes'}
    log("cells: " + ", ".join(f"{len(ids):,} {label[t]}" for t, _, ids in cells))
    if n_skipped:
        log(f"warning: {n_skipped:,} cells have an element type this script cannot "
            f"rebuild (polyhedra / hanging nodes) and were left out")

    # ---- write
    folder = out_base
    os.makedirs(folder, exist_ok=True)
    stem = os.path.basename(out_base)
    cell_zone = np.zeros(n_cells + 1, np.int32)
    for z, first, last in cell_zones:
        cell_zone[first:last + 1] = z
    kzone = cell_zone[np.flatnonzero(keep)]

    vol_items, used_names = [], set()

    def unique(name):
        base, i = safe(name), 2
        name = base
        while name in used_names:
            name, i = f"{base}_{i}", i + 1
        used_names.add(name)
        return name

    for z in np.unique(kzone):
        pieces = [(t, conn[kzone[ids] == z]) for t, conn, ids in cells]
        if not sum(len(c) for _, c in pieces):
            continue
        name = zname(int(z)) if z else 'cells'
        fn = unique('volume_' + name) + '.vtu'
        npts, ncl = write_vtu(os.path.join(folder, fn), xyz, pieces)
        vol_items.append((name, f"{stem}/{fn}"))
        log(f"wrote volume  {name}: {ncl:,} cells, {npts:,} points")

    surf_items = []
    for z in sorted(surf):
        pieces = [({3: VTK_TRI, 4: VTK_QUAD}[k], np.concatenate(v)) for k, v in sorted(surf[z].items())]
        name = zname(z)
        fn = unique('surface_' + name) + '.vtu'
        npts, ncl = write_vtu(os.path.join(folder, fn), xyz, pieces)
        surf_items.append((name, f"{stem}/{fn}"))
        log(f"wrote surface {name}: {ncl:,} faces")

    vtm = out_base + '.vtm'
    write_vtm(vtm, [('volume', vol_items), ('surfaces', surf_items)])
    log(f"finished. Open in ParaView: {vtm}")


if __name__ == '__main__':
    try:
        main()
    except MemoryError:
        log("out of memory")
        raise