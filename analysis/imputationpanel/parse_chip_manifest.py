# %%
import os
import re
from collections import Counter

workdir = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/chip"
path_manifest = f"{workdir}/InfiniumOmni2-5-8v1-5_A1.csv"
path_rsid_conversion = f"{workdir}/InfiniumOmni2-5-8v1-5_A1_b151_rsids.txt"

# %%
def parse_illumina_manifest(filepath, output_file, miss_file, return_chrpos=True, write=True):
    set_ilmn_id = set()
    set_chrpos = set()
    
    in_assay_section = False
    header = None
    col_index = {}

    snp_pattern = re.compile(r"\[([ACGT])\/([ACGT])\]")

    with open(filepath, mode='r') as f, open(output_file, mode='w', newline='') as fw, open(miss_file, mode='w', newline='') as fm:
        if write:
            fw.write(",".join(["IlmnID", "Name", "Chr", "Pos", "RefAllele", "AltAllele", "RefStrand"]) + "\n")

        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect start of assay section
            if line.startswith("[Assay]"):
                in_assay_section = True
                continue

            if not in_assay_section:
                continue

            # First line after [Assay] is the header
            if header is None:
                header = line.split(',')
                fm.write(",".join(header) + "\n")
                required_cols = ["IlmnID", "Name", "SNP", "Chr", "MapInfo", "RefStrand"]
                col_index = {col: header.index(col) for col in required_cols if col in header}
                continue

            # Process data lines
            parts = line.split(',')

            try:
                ilmn_id = parts[col_index["IlmnID"]].strip()
                set_ilmn_id.add(ilmn_id)
                ilmn_name = parts[col_index["Name"]].strip()
                snp_field = parts[col_index["SNP"]].strip()
                chr_val = parts[col_index["Chr"]].strip()
                pos_val = parts[col_index["MapInfo"]].strip()
                chr_pos_val  = chr_val + ":" + pos_val
                set_chrpos.add(chr_pos_val)
                ref_strand = parts[col_index["RefStrand"]].strip() if "RefStrand" in col_index else "-"

                # Parse SNP like [A/G]
                match = snp_pattern.match(snp_field)
                if match:
                    ref, alt = match.groups()
                else:
                    ref, alt = ("N", "N")

                chr_val = re.sub(r"[^0-9XYMTNAPar]+", "", chr_val.upper())

                if re.match(r"^(0|NA|XY|X|Y|MT|M|PAR)", chr_val):
                    continue
                
                if pos_val in ["0", "", "NA"]:
                    continue
                
                if write:
                    fw.write(",".join([ilmn_id, ilmn_name, chr_val, pos_val, ref, alt, ref_strand]) + "\n")

            except Exception as e:
                print(f"Problem line: {line[:80]} -> {e}")
                continue

    print(f"Finished parsing: {output_file}")
    
    if return_chrpos:
        return set_chrpos
    
    else:
        return set_ilmn_id

# %%
def get_ilmn_ids(path_manifest):
    set_ilmn_ids = set()
    header = None
    with open(path_manifest, mode="r") as fr:
        for line in fr:
            line = line.rstrip()
            if not line:
                continue
            
            if header is None:
                header = line.split(",")
                ilmn_idx = header.index("IlmnID")
                continue
            
            parts = line.split(",")
            ilmn_id = parts[ilmn_idx]
            set_ilmn_ids.add(ilmn_id)
    
    return set_ilmn_ids

# %%
path_outfile = os.path.join(workdir, "parsed_manifest.csv")
path_missfile = os.path.join(workdir, "missed_manifest.tsv")
set_original = parse_illumina_manifest(path_manifest, path_outfile, path_missfile, return_chrpos=False, write=False)
set_parsed = get_ilmn_ids(path_outfile)
set_miss = set_original - set_parsed
set_chrpos = parse_illumina_manifest(path_manifest, path_outfile, path_missfile, return_chrpos=True, write=True)

# %%
for i in range(1, 23, 1):
    path_taret_chrpos = f"{workdir}/InfiniumOmni2-5-8v1-5_A1.target.pos.chr{i}.list"
    with open(path_taret_chrpos, mode="w") as fw:
        for chrpos in set_chrpos:
            if str(chrpos).split(":")[0] == str(i):
                fw.write(f"chr{chrpos}\n")
            
# %%
