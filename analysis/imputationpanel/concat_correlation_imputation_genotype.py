# %%
#!/usr/bin/env python3
import glob
import gzip
import os

from joblib import Parallel, delayed


# %%
def concat_corr_files(dir_corr, chrom, panel):
    # Build pattern for input files
    pattern = f"chr{chrom}.dose_chunk_*.{panel}.corr.txt.gz"
    files = sorted(glob.glob(os.path.join(dir_corr, pattern)))

    if not files:
        raise FileNotFoundError(f"No matching files found for pattern: {pattern}")

    outfile = f"chr{chrom}.dose.{panel}.concat.corr.txt.gz"
    print(f"✅ Found {len(files)} files for {chrom} ({panel})")
    print(f"Output file: {outfile}")

    # Write header from first file
    with gzip.open(outfile, "wb") as out:
        with gzip.open(files[0], "rt") as f:
            header = f.readline()
            out.write(header.encode())

        # Append data lines from all files (skip header lines)
        for file in files:
            print(f"Appending: {file}")
            with gzip.open(file, "rt") as f:
                for line in f:
                    if line.startswith("CHROM"):
                        continue
                    out.write(line.encode())

    print(f"✅ Done. Merged {len(files)} chunk(s) into {outfile}")


dir_corr = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation/corr_results"
with Parallel(n_jobs=23) as parallel:
    parallel(delayed(concat_corr_files)(dir_corr, chrom=i, panel="4K") for i in range(1, 23, 1))