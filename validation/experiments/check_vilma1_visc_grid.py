#!/usr/bin/env python3
"""check_vilma1_visc_grid.py <densi.inp> <visc3d.nc>

Assert the contract VILMA1 does not check itself: a 3-D viscosity file must carry
exactly one level per viscoelastic element of the densi.inp it is run with, in
the same order.

VILMA1 maps the file's radius axis onto its elements BY INDEX, not by radius, and
silently ignores any levels past its element count. Measured (an earlier pass):
the 164-level Bagge-derived field run against the 129-element M3-L70-V01
densi.inp left VILMA1 reading levels 1..129 and stopping at r = 6193 km -- no
lithosphere, deep mantle stretched over the whole Earth, surface element 10^20.95
against the file's 10^29.99. It printed no warning.

Element geometry, as VILMA1 builds it: each densi.inp row after the first (the
fluid core, which carries no viscoelastic element) is `r_bot r_top dr rho mu`,
split into trunc((r_top-r_bot)/dr) elements identified by their centre radius.

ASSUMES liqcor = 1 and exactly one fluid-core row, which holds for every
densi.inp in this project: the core elements are stripped by gen_liqcor, so
CLIMBER-X's 222 elements become 164, the M3-L70-V01 file's 216 become 129, and
the VILMA2-PREM file's 215 become 128. With liqcor = 0, or a core split
over several rows, this count is wrong.
"""
import re
import subprocess
import sys


def fnum(s):
    """Fortran real: 3480.d03 / 1.0D+22 / plain decimal."""
    return float(s.replace("d", "e").replace("D", "E"))


def elements(densi):
    rows = [l.split() for l in open(densi).read().splitlines()[2:] if len(l.split()) >= 3]
    out = []
    for rb, rt, dr, *_ in rows[1:]:                     # skip the fluid core
        rb, rt, dr = fnum(rb), fnum(rt), fnum(dr)
        n = max(1, int((rt - rb) / dr))
        h = (rt - rb) / n
        out += [rb + (j + 0.5) * h for j in range(n)]
    return out


def main(densi, nc):
    els = elements(densi)
    # ncdump rather than a Python NetCDF module: one dependency fewer
    txt = subprocess.run(["ncdump", "-v", "radius", "-p", "9,17", nc], check=True,
                         capture_output=True, text=True).stdout
    data = txt[txt.index("data:"):]
    r = [float(x) for x in re.search(r"radius\s*=([^;]*);", data).group(1).split(",")]
    if len(els) != len(r):
        sys.exit(f"MISMATCH: {densi} has {len(els)} viscoelastic elements but {nc} has "
                 f"{len(r)} levels.\nVILMA1 maps this axis by index and would silently "
                 f"integrate the wrong Earth. Regenerate with "
                 f"experiments/make_vilma1_visc_ongrid.jl.")
    off = max(abs(a - b) for a, b in zip(els, r))
    if off > 1.0:
        sys.exit(f"MISMATCH: {nc} has the right number of levels but they are not "
                 f"{densi}'s element centres (max offset {off:.1f} m).")
    print(f"    viscosity grid check: {len(r)} levels = {len(els)} elements, "
          f"centres agree to {off:.3f} m")


if __name__ == "__main__":
    main(*sys.argv[1:3])
