#%%
import gzip


def get_sample_names(vcf_path):
    open_func = gzip.open if vcf_path.endswith(".bgz") else open
    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.strip().split('\t')
                return cols[9:]  # Samples start from 10th column

vcf_path = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/grf_vcf/VCF_ChipOnly_MatchedOnly/chr21.omniChipOnly.vcf"
sample_path = "/BiO/Access/kyungwhan1998/genome/shapeit/Resources/Data/grf_vcf/VCF_ChipOnly_MatchedOnly/samples.list"

samples = get_sample_names(vcf_path)


with open(sample_path, mode="w") as fw:
    for sample in samples:
        fw.write(sample+"\n")