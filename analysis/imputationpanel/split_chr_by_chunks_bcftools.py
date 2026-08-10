#!/usr/bin/env python3

import os
import subprocess
import sys

import pandas as pd
from joblib import Parallel, delayed


def run_cmd(cmd):
    """Run a shell command safely."""
    print(f"[CMD] {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)

def extract_region(in_vcf, chunk_id, region, outdir):
    """Extract VCF subset for a specific region using bcftools."""
    outprefix = os.path.basename(in_vcf).replace(".vcf.gz", "")
    out_vcf = os.path.join(outdir, f"{outprefix}_chunk_{chunk_id}.vcf.gz")
    cmd = f"bcftools view -r {region} -Oz -o {out_vcf} {in_vcf} && tabix -p vcf {out_vcf}"
    run_cmd(cmd)

    return out_vcf

def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <input.vcf.gz> <regions.tsv> <output_dir> <threads>")
        sys.exit(1)

    vcf = sys.argv[1]
    region_file = sys.argv[2]
    outdir = sys.argv[3]
    threads = sys.argv[4]

    os.makedirs(outdir, exist_ok=True)

    # Read region table (TSV format)
    df = pd.read_csv(region_file, sep=r"\s+")

    # Parallel execution
    print(f"Starting joblib parallel with {threads} workers...", flush=True)
    results = Parallel(n_jobs=int(threads))(
        delayed(extract_region)(vcf, row["chunk"], row["input_region"], outdir)
        for _, row in df.iterrows()
    )

    print("✅ Done: all chunks generated and indexed.")

if __name__ == "__main__":
    main()
