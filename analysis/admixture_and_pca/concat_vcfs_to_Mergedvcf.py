# %%
import glob
import os
import subprocess


# %%
def create_concat_file_list(path_file_list, dir_vcf, pattern=".recalibrated_variants.vcf.gz", list_remove_tag=["chrX", "chrY", "others"]):
    file_list = glob.glob(f"{dir_vcf}/*{pattern}")
    filter_file_list = list(filter(lambda x: os.path.basename(x).split("_")[-2].split(".")[0] not in list_remove_tag, file_list))
    with open(path_file_list, mode="w") as fw:
        for file in filter_file_list:
            fw.write(file+"\n")
    
# %%
bcftools = "/BiO/Research/Korea10KGenome/Resources/Tools/bcftools-1.20/bcftools"
dir_vcf = "/BiO/Access/kyungwhan1998/Korea10KGenome/Resources/External_Genome_Data/1KGP/1KGP_30x_GRCh38"
path_file_list = f"{dir_vcf}/file_vcf_list.txt"
create_concat_file_list(path_file_list, dir_vcf)
outfile = f"{dir_vcf}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_merged_chr.recalibrated_variants.annotated.vcf.gz"

cmd = f"{bcftools} concat --file-list {path_file_list} -Oz -o {outfile}"

# %%
subprocess.run(cmd, shell=True)