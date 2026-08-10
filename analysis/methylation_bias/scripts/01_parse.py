#!/usr/bin/env python3
"""
Parse per_cpg_beta_by_genotype.with_values.filtered.tsv into a compact binary store.

Output: data/betas.npz
    chrom, pos            : per-site identifiers
    n0, n1, n2            : genotype counts
    off0, off1, off2      : start offsets into the flat value arrays
    v0, v1, v2            : flat float32 arrays of beta values, concatenated per site
    af                    : alternate allele frequency = (n1 + 2*n2) / (2*N)
"""
import argparse
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "per_cpg_beta_by_genotype.with_values.filtered.tsv"
OUT = ROOT / "data" / "betas.npz"


def main(src=SRC):
    chrom, pos = [], []
    n0l, n1l, n2l = [], [], []
    v0, v1, v2 = [], [], []

    with open(src) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            chrom.append(p[0])
            pos.append(int(p[1]))
            n0l.append(int(p[2]))
            n1l.append(int(p[6]))
            n2l.append(int(p[10]))
            for col, acc in ((14, v0), (15, v1), (16, v2)):
                s = p[col]
                if s:
                    acc.append(np.fromstring(s, sep=",", dtype=np.float64).astype(np.float32))

    n0 = np.array(n0l, dtype=np.int32)
    n1 = np.array(n1l, dtype=np.int32)
    n2 = np.array(n2l, dtype=np.int32)
    N = n0.astype(np.int64) + n1 + n2
    af = (n1 + 2.0 * n2) / (2.0 * N)

    def offsets(n):
        o = np.zeros(len(n) + 1, dtype=np.int64)
        np.cumsum(n, out=o[1:])
        return o

    np.savez_compressed(
        OUT,
        chrom=np.array(chrom),
        pos=np.array(pos, dtype=np.int64),
        n0=n0, n1=n1, n2=n2, af=af,
        off0=offsets(n0), off1=offsets(n1), off2=offsets(n2),
        v0=np.concatenate(v0), v1=np.concatenate(v1), v2=np.concatenate(v2),
    )
    print(f"wrote {OUT}  sites={len(n0)}  N_median={np.median(N):.0f}", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default=str(SRC),
                    help="input TSV (default: the study filename in the repo root)")
    main(ap.parse_args().input)
