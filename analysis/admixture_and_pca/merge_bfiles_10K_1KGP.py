# %%
import glob
import os
import subprocess

# %%
plink2 = "/BiO/Share/Tool/plink2"
indir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data"

# %%
def run_write_snplist(plink2, bfile_in, snpfile_out):
    cmd = f"{plink2} --bfile {bfile_in} --write-snplist --out {snpfile_out}"
    
    return cmd

# %%
def run_extract_intersect(plink2, bfile_in, snpfile_in, bfile_out):
    cmd = f"{plink2} --bfile {bfile_in} \
                --extract-intersect {snpfile_in}\
                --make-bed \
                --out {bfile_out}"
    return cmd

# %%
def run_bmerge(plink2, bfile_in1, bfile_in2, bfile_out, test=False):
    if test:
        cmd = f"{plink2} --bfile {bfile_in1} --bmerge {bfile_in2} --merge-mode 6 --out {bfile_out}.missnp"
    else:
        cmd = f"{plink2} --bfile {bfile_in1} --bmerge {bfile_in2} --make-bed --out {bfile_out}"
    
    return cmd

# %%
def run_filp(plink1, bfile_in2, missnp, bfile_out):
    cmd = f"{plink1} --bfile {bfile_in2} --flip {missnp} --make-bed --out {bfile_out}"
    
    return cmd

# %%
def run_recode_vcf(plink1, bfile_in, vcf_out):
    cmd = f"{plink1} --bfile {bfile_in} --recode vcf --out {vcf_out}"
    
    return cmd

# %%
def run_exclude(plink1, bfile_in, missnp, bfile_out):
    cmd = f"{plink1} --bfile {bfile_in} --exclude {missnp} --make-bed --out {bfile_out}"
    
    return cmd

# %%
indir_10K = os.path.join(indir, "ku10k")
indir_1KGP = os.path.join(indir, "1kgp")
outdir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_10K = list(map(lambda x: ".".join(x.split(".")[:-1]), glob.glob(f"{indir_10K}/*.bed")))
bfile_10K = list(filter(lambda x: "extracted" not in x, bfile_10K))[0]
bfile_1KGP = list(map(lambda x: ".".join(x.split(".")[:-1]), glob.glob(f"{indir_1KGP}/*.bed")))
bfile_1KGP = list(filter(lambda x: "extracted" not in x, bfile_1KGP))[0]
os.makedirs(outdir, exist_ok=True)

snpfile_10K = os.path.join(outdir, "10K")
cmd1 = run_write_snplist(plink2, bfile_10K, snpfile_10K)
print(cmd1)
snpfile_1KGP = os.path.join(outdir, "1KGP")
cmd2 = run_write_snplist(plink2, bfile_1KGP, snpfile_1KGP)
print(cmd2)

# %%
bfile_10K_1KGP_snp_extract = bfile_10K + ".1kgp_snp_extracted"
cmd1 = run_extract_intersect(plink2, bfile_10K, snpfile_1KGP+".snplist", bfile_10K_1KGP_snp_extract)
print(cmd1)
bfile_1KGP_10K_snp_extract = bfile_1KGP + ".ku10k_snp_extracted"
cmd2 = run_extract_intersect(plink2, bfile_1KGP, snpfile_10K+".snplist", bfile_1KGP_10K_snp_extract)
print(cmd2)

# %%
plink1 = "/BiO/Share/Tool/plink"
bfile_10K_1KGP_snp_extract = bfile_10K + ".1kgp_snp_extracted"
bfile_1KGP_10K_snp_extract = bfile_1KGP + ".ku10k_snp_extracted"
bfile_out = os.path.join(outdir, "Merged_ku10k_1kgp.extract_overlap_snps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
cmd = run_bmerge(plink1, bfile_10K_1KGP_snp_extract, bfile_1KGP_10K_snp_extract, bfile_out, test=True)
print(cmd)

# %%
missnp = f"{bfile_out}.missnp.missnp"
bfile_1KGP_10K_snp_extract_flipped = bfile_1KGP_10K_snp_extract + ".flipped"
cmd = run_filp(plink1, bfile_1KGP_10K_snp_extract, missnp, bfile_1KGP_10K_snp_extract_flipped)
print(cmd)

# %%
bfile_out = os.path.join(outdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
cmd = run_bmerge(plink1, bfile_10K_1KGP_snp_extract, bfile_1KGP_10K_snp_extract_flipped, bfile_out, test=True)
print(cmd)

# %%
vcf_out = f"{bfile_10K_1KGP_snp_extract}.vcf"
cmd1 = run_recode_vcf(plink1, bfile_10K_1KGP_snp_extract, vcf_out)
print(cmd1)
vcf_out = f"{bfile_1KGP_10K_snp_extract_flipped}.vcf"
cmd2 = run_recode_vcf(plink1, bfile_1KGP_10K_snp_extract_flipped, vcf_out)
print(cmd2)

# %%
bfile_10K_1KGP_snp_extract_missnp_exclude = f"{bfile_10K_1KGP_snp_extract}.missnp_excluded"
bfile_1KGP_10K_snp_extract_flipped_missnp_exclude = f"{bfile_1KGP_10K_snp_extract}.missnp_excluded"
cmd1 = run_exclude(plink1, bfile_10K_1KGP_snp_extract, missnp, bfile_10K_1KGP_snp_extract_missnp_exclude)
cmd2 = run_exclude(plink1, bfile_1KGP_10K_snp_extract_flipped, missnp, bfile_1KGP_10K_snp_extract_flipped_missnp_exclude)
print(cmd1)
print(cmd2)

# %%
bfile_out = os.path.join(outdir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
cmd = run_bmerge(plink1, bfile_10K_1KGP_snp_extract_missnp_exclude, bfile_1KGP_10K_snp_extract_flipped_missnp_exclude, bfile_out, test=False)
print(cmd)

# %%
plink2 = "/BiO/Share/Tool/plink2"
indir = "/BiO/Access/kyungwhan1998/genome/admixture/Resources/Data/ku10k_1KGP"
bfile_in = list(map(lambda x: ".".join(x.split(".")[:-1]), glob.glob(f"{indir}/*.bed")))[0]
bfile_out = os.path.join(indir, "Merged_ku10k_1kgp.extract_overlap_snps.flipped_nonoverlap_snps.excluded_missnps.postmerge_QC_filtered.Merged_chr.biallelic.Autosome.varname.geno_0.01.mind_0.1.hwe_1e6.het_3std.kinship_3rd")
cmd = f"{plink2} --bfile {bfile_in} --geno 0.01 --maf 0.01 --hwe 1e-6 --make-bed --out {bfile_out}"
print(cmd)
