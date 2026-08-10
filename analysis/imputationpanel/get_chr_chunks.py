# %%
import os

from BioParser import fileparser


# %%
def get_dict_contig_length(file_vcf):
    vcf = fileparser.VCF(file_vcf, is_gzip=True, metadata=True)
    metadata_vcf = vcf.metadata
    contigs_vcf = list(filter(lambda x: str(x).startswith("##contig"), metadata_vcf))
    dict_contig_len = dict()
    for contig_vcf in contigs_vcf:
        contig_content = contig_vcf.split(",")
        contig_id = contig_content[0].split("=")[-1]
        contig_length = contig_content[1].split("=")[-1].replace(">", "")
        
        if len(contig_id.split("_")) > 1:
            break
        
        else:
            dict_contig_len[contig_id] = int(contig_length)
        
    return dict_contig_len

# %%
def write_chunk_coordinates(file_chunk, dict_contig_len, chunk_interval, chunk_overlap):
    chunk_id = 0
    with open(file_chunk, mode="w") as fw:
        for contig_id, contig_length in dict_contig_len.items():
            chunk_start = 1
            while chunk_start <= contig_length:
                chunk_end = min(chunk_start + chunk_interval - 1, contig_length)
                if chunk_end > contig_length:
                    chunk_end = contig_length
                fw.write("\t".join([str(chunk_id), str(contig_id), str(chunk_start), str(chunk_end)]) + "\n")
                if chunk_end == contig_length:
                    break
                chunk_start += (chunk_interval - chunk_overlap) 
                chunk_id += 1
                
# %%
file_vcf = "/BiO/Access/kyungwhan1998/genome/shapeit/chr22.biallelic.vcf.gz"
file_chunk = "/BiO/Access/kyungwhan1998/genome/shapeit/chunks.coordinates.txt"
Mbp = 1000000
chunk_interval = (6*Mbp)
chunk_overlap = (2*Mbp)
dict_contig_len = get_dict_contig_length(file_vcf)
write_chunk_coordinates(file_chunk, dict_contig_len, chunk_interval, chunk_overlap)