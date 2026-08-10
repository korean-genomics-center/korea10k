#!/usr/bin/env python3
"""
Generate a synthetic input file in the same format as the real one, so that the
pipeline can be run end to end without access to individual-level cohort data.

THIS IS A FORMAT DEMONSTRATION, NOT THE STUDY DATA. It reproduces the schema and
the broad shape of the real panel — an allele-frequency spectrum skewed to rare
variants, beta values that fall with alternate allele dosage, hyper- and
hypomethylated sites, and the 0/1 boundary pile-up — but the values are drawn
from a parametric model and carry no biological information. Numbers obtained
from it will not match the published figure; only the code path is the same.

To reproduce the published figure, obtain the real dataset (see data/README.md)
and skip this script.

Usage:
    python3 scripts/make_synthetic_input.py --sites 2000 --samples 1208
"""
import argparse
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT = ROOT / "per_cpg_beta_by_genotype.with_values.filtered.tsv"

HEADER = ("chr\tpos\tn_gt0\tmean_gt0\tvar_gt0\tmed_gt0\t"
          "n_gt1\tmean_gt1\tvar_gt1\tmed_gt1\t"
          "n_gt2\tmean_gt2\tvar_gt2\tmed_gt2\tbeta_gt0\tbeta_gt1\tbeta_gt2\n")


def beta_draw(rng, mean, conc, size):
    """Beta-distributed values around `mean`, clipped to [0, 1] and rounded
    to four decimals, matching the precision of the real file."""
    mean = np.clip(mean, 1e-3, 1 - 1e-3)
    a, b = mean * conc, (1 - mean) * conc
    return np.round(np.clip(rng.beta(a, b, size), 0.0, 1.0), 4)


def main(args):
    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    n_sites, N = args.sites, args.samples

    # allele frequencies skewed toward rare, as in a real variant panel
    af = np.clip(rng.beta(0.7, 2.2, n_sites), 0.001, 0.995)
    # a majority of hypermethylated sites, a minority hypomethylated
    hyper = rng.random(n_sites) < 0.75
    base = np.where(hyper, rng.uniform(0.80, 0.98, n_sites),
                    rng.uniform(0.02, 0.20, n_sites))

    with open(out, "w") as fh:
        fh.write(HEADER)
        for i in range(n_sites):
            p = af[i]
            # genotype counts from HWE with multinomial sampling noise
            n0, n1, n2 = rng.multinomial(N, [(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
            if min(n0, n1, n2) < 1:          # the real panel requires all three
                n0, n1, n2 = max(n0, 1), max(n1, 1), max(n2, 1)

            # methylation falls toward the opposite extreme with alt dosage
            far = 0.01 if hyper[i] else 0.99
            m0 = base[i]
            m2 = far + (rng.normal(0, 0.02) if far < 0.5 else rng.normal(0, 0.02))
            m2 = float(np.clip(m2, 0.0, 1.0))
            # heterozygotes near the midpoint, with a small dominance deviation
            m1 = float(np.clip((m0 + m2) / 2 + rng.normal(0.01, 0.03), 0.0, 1.0))

            v0 = beta_draw(rng, m0, rng.uniform(60, 220), n0)
            v1 = beta_draw(rng, m1, rng.uniform(25, 90), n1)
            v2 = beta_draw(rng, m2, rng.uniform(60, 220), n2)

            cols = [f"chr{rng.integers(1, 23)}", str(int(rng.integers(1e5, 2.5e8)))]
            for v in (v0, v1, v2):
                cols += [str(len(v)), f"{v.mean():.5f}", f"{v.var():.5f}",
                         f"{np.median(v):.5f}"]
            cols += [",".join(f"{x:.4f}" for x in v) for v in (v0, v1, v2)]
            fh.write("\t".join(cols) + "\n")

    mb = out.stat().st_size / 1e6
    print(f"wrote {out}  ({n_sites} synthetic sites, {N} samples, {mb:.1f} MB)")
    print("NOTE: synthetic data for pipeline testing only; results will not "
          "match the published figure.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", type=int, default=2000)
    ap.add_argument("--samples", type=int, default=1208)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    main(ap.parse_args())
