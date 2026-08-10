# %%
import gzip
import os

from joblib import Parallel, delayed


# %%
def inner_merge_by_variant(i, panel, R_template, maf_template, R_maf_template, header,
                           var_col_maf=1, freq_col_mac=-1, freq_col_maf=-2, var_col_r=1, r_col_r=-1):
    """
    Merge R and MAF files by matching variants.
    """
    R_chr = R_template.format(i=i, panel=panel)
    maf_chr = maf_template.format(i=i)
    R_maf_chr = R_maf_template.format(i=i, panel=panel)

    print(f"Merging chr{i} ...")

    # Read MAF into dictionaries
    maf_dict = {}
    mac_dict = {}
    with open(maf_chr, "r") as f_maf:
        f_maf.readline()
        for line in f_maf:
            cols = line.rstrip("\n").split("\t")
            variant = int(cols[var_col_maf].split("_")[-1])
            maf_dict[variant] = cols[freq_col_maf]
            mac_dict[variant] = cols[freq_col_mac]

    # Merge R with MAF (gzipped)
    with gzip.open(R_chr, "rt") as f_R, gzip.open(R_maf_chr, "wt") as fw:
        fw.write("\t".join(header) + "\n")
        f_R.readline()

        for line in f_R:
            cols = line.rstrip("\n").split("\t")
            variant = int(cols[var_col_r])
            r = cols[r_col_r]
            if variant in maf_dict:
                fw.write("\t".join([str(variant), r, mac_dict[variant], maf_dict[variant]]) + "\n")

    print(f"✅ chr{i} merged → {R_maf_chr}")

# %%
def categorize_variant(ac, af):
    """Return category name based on AC and AF thresholds."""
    try:
        ac = int(ac)
        af = float(af)
    except ValueError:
        return "NA"

    if ac == 1:
        return "Singleton"
    elif ac == 2:
        return "Doubleton"
    elif 2 < ac and af <= 0.0005:
        return "Extremely Rare"
    elif 0.0005 < af <= 0.001:
        return "Very Rare"
    elif 0.001 < af <= 0.01:
        return "Rare"
    elif 0.01 < af <= 0.05:
        return "Common"
    elif af > 0.05:
        return "Very Common"
    else:
        return "NA"

# %%
def add_variant_category(i, panel, R_maf_template, R_maf_cat_template, header):
    """
    Add AF category column to the merged R + MAF file.
    """
    R_maf_chr = R_maf_template.format(i=i, panel=panel)
    R_maf_cat_chr = R_maf_cat_template.format(i=i, panel=panel)

    print(f"Processing chr{i}...")

    with gzip.open(R_maf_chr, "rt") as fin, gzip.open(R_maf_cat_chr, "wt") as fout:
        header_line = fin.readline().rstrip("\n")
        # Write new header with AF category
        fout.write("\t".join(header + ["ALT_CNT", "AF_Category"]) + "\n")

        for line in fin:
            cols = line.rstrip("\n").split("\t")
            af = float(cols[3])
            ac = int(float(cols[2]) * float(cols[3]))
            category = categorize_variant(ac, af)
            fout.write(line.rstrip("\n") + "\t" + str(ac) + "\t" + category + "\n")

    print(f"✅ chr{i} categorized → {R_maf_cat_chr}")

# %%
# Paths and templates
workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Results/imputation"

R_template = os.path.join(workdir, "corr_results/chr{i}.dose.{panel}.concat.corr.txt.gz")
maf_template = os.path.join(workdir, "maf_info/chr{i}_maf_10K.txt")
R_maf_template = os.path.join(workdir, "corr_results/chr{i}.dose.{panel}.concat.corr.maf_10K_added.txt.gz")
R_maf_cat_template = os.path.join(workdir, "corr_results/chr{i}.dose.{panel}.concat.corr.maf_10K_cat_added.txt.gz")

header = ["Variant", "R", "OBS_CNT", "ALT_FREQ"]
panel = "10K"

# %%
def run(i):
    inner_merge_by_variant(i, panel, R_template, maf_template, R_maf_template, header)
    add_variant_category(i, panel, R_maf_template, R_maf_cat_template, header)

#
with Parallel(n_jobs=23) as parallel:
    parallel(delayed(run)(i) for i in range(1, 23))
