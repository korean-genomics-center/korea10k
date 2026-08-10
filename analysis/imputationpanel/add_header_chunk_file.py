# %%
for num in range(1, 23, 1):
    file_chunk_coord = f"/BiO/Share/Tool/shapeit5/resources/chunks/b38/20cM/chunks_chr{num}.txt"
    file_out = f"/BiO/Share/Tool/shapeit5/resources/chunks/b38/20cM/chunks_chr{num}.header.incl.txt"
    with open(file_chunk_coord, mode='r') as fr, open(file_out, mode="w") as fw:
        fw.write("\t".join(["chunk", "chr", "scaffold_region", "input_region"]) + "\n")
        for line in fr:
            record = line.rstrip().split("\t")[:4]
            record_rmv_chr = list(map(lambda x: x.replace("chr", ""), record))
            fw.write("\t".join(record_rmv_chr)+"\n")
            
# %%
file_merged = "/BiO/Share/Tool/shapeit5/resources/chunks/b38/20cM/chunks.header.incl.txt"

with open(file_merged, mode="w") as fw:
    fw.write("\t".join(["chunk", "chr", "scaffold_region", "input_region"]) + "\n")
    for num in range(1, 23, 1):
        file_chunks = f"/BiO/Share/Tool/shapeit5/resources/chunks/b38/20cM/chunks_chr{num}.header.incl.txt"
        with open(file_chunks, mode='r') as fr:
            for line in fr:
                if str(line).startswith("chunk"):
                    continue
                else:
                    record = line.rstrip().split("\t")[:4]
                    fw.write("\t".join(record)+"\n")
# %%
