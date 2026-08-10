# %%
import re

import pandas as pd

# %%
# -----------------------------
# 1. Load input HLA typing table
# -----------------------------
hla_path = "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/20140702_hla_diversity.txt"
df = pd.read_csv(hla_path, delim_whitespace=True)

# %%
# -----------------------------
# 2. Define loci and output columns
# -----------------------------
loci = ["A", "B", "C", "DRB1", "DQB1"]
output_cols = ["SampleID"] + [f"{locus}_{i}" for locus in loci for i in [1, 2]]

# %%
# -----------------------------
# 3. Helper functions
# -----------------------------
def normalize_allele_2digit(allele, locus):
    """Add locus prefix and truncate to 2-digit resolution (xx:xx)."""
    allele = str(allele).strip()
    if not allele or allele.upper() == "NA":
        return ""
    # keep only first two fields after splitting ':'
    parts = allele.split(":")[:2]
    allele_2digit = ":".join(parts)
    if not allele_2digit.startswith(f"{locus}*"):
        allele_2digit = f"{locus}*{allele_2digit}"
    return allele_2digit

def collapse_alleles_2digit(value, locus):
    """Collapse ambiguous alleles to 2-digit forms, keep all unique if different."""
    if pd.isna(value) or value == "":
        return [""]
    parts = [normalize_allele_2digit(v, locus) for v in re.split(r"[ /]", str(value)) if v.strip()]
    # remove duplicates
    return sorted(set(parts))

# %%
# -----------------------------
# 4. Parse and standardize
# -----------------------------
output_rows = []
for _, row in df.iterrows():
    record = {"SampleID": row["id"]}

    for locus in loci:
        a1 = collapse_alleles_2digit(row[locus], locus)
        a2 = collapse_alleles_2digit(row[f"{locus}.1"], locus)

        record[f"{locus}_1"] = "/".join(a1) if a1 else ""
        record[f"{locus}_2"] = "/".join(a2) if a2 else ""

    output_rows.append(record)

hla_clean = pd.DataFrame(output_rows)[output_cols]

# %%
# -----------------------------
# 5. Remove samples with missing alleles (*00:00)
# -----------------------------
mask_valid = ~hla_clean.apply(lambda row: any("*0000" in str(x) for x in row[1:]), axis=1)
hla_clean_filtered = hla_clean[mask_valid].copy()

# %%
# -----------------------------
# 6. Save filtered table to TSV
# -----------------------------
hla_clean_filtered.to_csv(
    "/BiO/Access/kyungwhan1998/genome/hla/Resources/Data/1KGP_1113samples_HLA.tsv",
    sep="\t",
    index=False
)

print(f"Filtered samples: {hla_clean_filtered.shape[0]} / {hla_clean.shape[0]}")
print(hla_clean_filtered.head())
