#!/usr/bin/env python3
"""
Pick the sites the simulation runs on, stratified by allele frequency.

AF is the axis that matters here because a single cohort of ~1,208 people has to
supply both the gt0 pool (which sets the reachable cohort size n) and the variant
carriers used as donors. Those two shrink in opposite directions as AF moves, so
any statement about bias has to be read within an AF stratum, never pooled.

Sites are required to have at least 200 gt0 samples (the smallest cohort in the
grid) and at least one carrier of each donor genotype, then sampled evenly
across AF bins.
"""
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BINS = np.array([0.0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 1.0])
PER_BIN = 220
SEED = 11


def main():
    d = np.load(ROOT / "data" / "betas.npz", allow_pickle=False)
    n0, n1, n2, af = d["n0"], d["n1"], d["n2"], d["af"]

    eligible = (n0 >= 200) & (n1 >= 1) & (n2 >= 1)
    rng = np.random.default_rng(SEED)
    binid = np.digitize(af, BINS[1:-1], right=True)

    chosen = []
    print(f"{'AF bin':>14}  {'eligible':>8}  {'chosen':>6}  {'med n_gt0':>9}  {'med n_gt2':>9}")
    for b in range(len(BINS) - 1):
        cand = np.where(eligible & (binid == b))[0]
        take = cand if len(cand) <= PER_BIN else rng.choice(cand, PER_BIN, replace=False)
        chosen.append(take)
        lab = f"({BINS[b]:.2f},{BINS[b+1]:.2f}]"
        print(f"{lab:>14}  {len(cand):>8}  {len(take):>6}  "
              f"{np.median(n0[take]):>9.0f}  {np.median(n2[take]):>9.0f}")

    sel = np.sort(np.concatenate(chosen))
    np.save(ROOT / "data" / "site_selection.npy", sel)
    print(f"\ntotal sites selected: {len(sel)}")


if __name__ == "__main__":
    main()
