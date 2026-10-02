"""
Parsing of (possibly multi-geometry) XYZ files.

The XYZ format is positional, not content-addressable:

    <natoms>
    <comment>              -- ARBITRARY: may be blank, may contain spaces,
    <symbol> <x> <y> <z>      may look exactly like an atom line
    ...                       (natoms of these)
    <natoms>               -- next frame, and so on

The comment line therefore has to be consumed BY POSITION.  Filtering lines
by what they look like -- e.g. dropping blank and single-token lines and
keeping the rest -- silently works for comments like "4.5" or
"Hexatriene,^1A_g,CC3" and then fails the moment a comment contains a space,
at which point it is parsed as an atom and the run either crashes in
float() or, worse, reshapes into a wrong geometry.
"""

import numpy as np


class XYZError(Exception):
    """Malformed XYZ input, reported with file and line number."""


def _err(path, lineno, msg):
    raise XYZError(f'{path}, line {lineno}: {msg}')


def read_frames(path, require_uniform=True):
    """Parse an XYZ file into a list of (symbols, coords) frames.

    symbols is a list of str; coords is an (natoms, 3) float array.  Blank
    lines between frames, a missing trailing newline, and comment lines of
    any content (including empty) are all tolerated.  Columns beyond the
    first four on an atom line are ignored, which some writers append.

    require_uniform: every frame must carry the same atoms in the same
    order, as a scan or trajectory must.
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    frames = []
    i, n = 0, len(lines)
    while True:
        # between frames, blank lines carry no meaning -- skip to the count
        while i < n and lines[i].strip() == '':
            i += 1
        if i >= n:
            break

        count_lineno = i + 1
        tok = lines[i].split()
        if len(tok) != 1:
            _err(path, count_lineno,
                 f'expected an atom count, found {len(tok)} fields: '
                 f'{lines[i].strip()!r}')
        try:
            natm = int(tok[0])
        except ValueError:
            _err(path, count_lineno,
                 f'expected an atom count, found {tok[0]!r}')
        if natm <= 0:
            _err(path, count_lineno, f'atom count must be positive, got {natm}')
        i += 1

        # the comment line is consumed positionally and never interpreted
        if i >= n:
            _err(path, count_lineno + 1,
                 'file ends where a comment line was expected')
        i += 1

        if i + natm > n:
            _err(path, i + 1,
                 f'frame starting at line {count_lineno} declares {natm} atoms '
                 f'but only {n - i} lines remain')

        syms, crds = [], []
        for k in range(natm):
            lineno = i + k + 1
            fld = lines[i + k].split()
            if len(fld) < 4:
                _err(path, lineno,
                     f'expected "symbol x y z", found {len(fld)} fields: '
                     f'{lines[i + k].strip()!r}')
            try:
                crds.append([float(x) for x in fld[1:4]])
            except ValueError:
                _err(path, lineno,
                     f'could not read coordinates from {lines[i+k].strip()!r}')
            syms.append(fld[0])
        i += natm
        frames.append((syms, np.array(crds, dtype=float)))

    if not frames:
        raise XYZError(f'{path}: no geometries found')

    if require_uniform:
        ref = frames[0][0]
        for j, (syms, _) in enumerate(frames[1:], start=2):
            if len(syms) != len(ref):
                raise XYZError(f'{path}: frame {j} has {len(syms)} atoms, '
                               f'frame 1 has {len(ref)}')
            bad = [(a, b) for a, b in zip(syms, ref) if a.lower() != b.lower()]
            if bad:
                raise XYZError(f'{path}: frame {j} atom ordering differs from '
                               f'frame 1 (e.g. {bad[0][0]!r} vs {bad[0][1]!r})')
    return frames
