# korea10k

Korea10K: 10,239 whole genomes with multiomic and clinical health information as
the Korean multiomics reference dataset (An et al., 2026, In submission).

This repository holds the analysis and figure-generation code for the study. It
is a code release, not a runnable end-to-end pipeline: the scripts expect the
project's genotype, methylation and metadata files to already exist on disk, and
they invoke external tools (PLINK 2, ADMIXTURE, SHAPEIT5, Minimac3/4, bcftools,
GATK, GLIMPSE2) that must be installed separately.

---

## Repository layout

```
korea10k/
├── korea10k/            # shared helper package (path configuration)
├── analysis/            # analysis pipelines, grouped by topic
│   ├── admixture_and_pca/
│   ├── depth/
│   ├── hla/
│   ├── imputationpanel/
│   ├── methylation_bias/
│   └── variantdiscovery/
├── paper/               # scripts that produce the manuscript figures/tables
│   ├── Supplementary_Figures/
│   └── Supplementary_Tables/
├── .env.example         # template for the site-specific path configuration
└── pyproject.toml
```

---

## Setup

### 1. Install the helper package

The scripts import a shared configuration module, so install the repository once
in your environment (editable, so edits are picked up immediately):

```bash
pip install -e .
```

### 2. Configure the data locations

**No absolute path to any data, tool or reference is stored in this
repository.** Every such location is read from an environment variable at run
time. Copy the template and fill in the directories for your machine:

```bash
cp .env.example .env
$EDITOR .env
```

`.env` is git-ignored, so real paths never enter version control. You can also
export the variables in your shell instead — an exported value always takes
precedence over the file.

| Variable | Points at |
| --- | --- |
| `KOREA10K_WORK_DIR` | Analysis workspace holding the per-analysis working directories (`admixture/`, `shapeit/`, `hla/`, `depthCoverage/`, `pca/`, `variant/`, `paper/`, …) |
| `KOREA10K_PROJECT_DIR` | Korea10K project root holding `Resources/`, `Results/` and `Analysis/` |
| `KOREA10K_TOOL_DIR` | Shared tool installations (PLINK 2, ADMIXTURE, SHAPEIT5, Minimac3/4, GLIMPSE2, GATK, …) |
| `KOREA10K_REF_DIR` | Reference genomes and liftover chain files (`hg38.fa`, `hg19.fa`, `hg19ToHg38.over.chain.gz`, …) |
| `KOREA10K_STORE_DIR` | Bulk storage for raw sequencing and per-sample intermediate data |
| `KOREA10K_CONDA_ENV_DIR` | Conda environment root providing the binaries the pipeline shells out to (`tabix`, `bcftools`, …) |
| `KOREF_DIR` | KOREF personal multiomics reference project root |

Only the roots used by the scripts you actually run need to be set. Importing a
root that is unset raises `korea10k.config.MissingPathError` with a message
naming the variable, rather than silently resolving to a wrong location.

### 3. Use it from a script

```python
from korea10k.config import PROJECT_DIR, TOOL_DIR

plink = f"{TOOL_DIR}/plink2"
psam  = f"{PROJECT_DIR}/Results/.../chr21.biallelic.psam"
```

---

## Running the scripts

Most scripts are written as VS Code / Jupyter *interactive* scripts: they are
divided into `# %%` cells, run top to bottom, and take their parameters from
module-level constants near the top of the file rather than from the command
line. To reproduce a step, open the file, adjust the constants (chromosome,
sample list, panel size, thresholds) and run it.

Scripts whose name begins with `run_`, `make_` or `calc_` generally shell out to
an external tool, and several of them submit the work to a Grid Engine cluster
via `qsub`; adapt the submission helpers if your scheduler differs.

---

## `analysis/`

### `admixture_and_pca/` — population structure

Merges Korea10K with 1KGP, applies sample/variant QC, and produces the PCA and
ADMIXTURE results.

| Stage | Scripts |
| --- | --- |
| Prepare genotypes | `concat_vcfs_to_Mergedvcf.py`, `convert_pgen_to_bed.py`, `merge_bfiles_10K_1KGP.py` |
| QC and pruning | `filter_samples_plink.py`, `filter_variants_plink.py`, `prune_LDs.py`, `get_samples_list.py` |
| PCA | `calc_pca.py`, `calc_pca_for_projection.py`, `calc_k_clusters_plink.py` |
| PCA plots | `draw_pca.py`, `draw_pca_postQC.py`, `draw_pca_cluster_annot.py`, `draw_pca_projection.py`, `draw_cpca.py` |
| ADMIXTURE | `run_admixture.py`, `plot_cv_error_admixture_choose_K.py`, `prepare_input_admixture_plot.py`, `draw_admixture_plot.py` |

`plot_cv_error_admixture_choose_K.py` plots the cross-validation error across
values of *K* and is what the choice of *K* is based on.

### `depth/` — sequencing depth and cohort metadata

- `calculate_depth.py` — converts per-sample mapped-base counts into average
  genome-wide depth and joins the sequencing platform from the metadata table.
- `plot_piechart_metadata.py` — cohort composition summary plots.
- `process_data_10239.py` — assembles the final 10,239-sample metadata table.

### `hla/` — HLA typing comparison

- `convert_1kgp_ku10k_format_hla.py` — puts the 1KGP and Korea10K HLA calls into
  a common format.
- `merge_1kgp_ku10k_hla_meta.py` — joins the calls with population metadata.
- `get_pop_hla_subtypes.py` — per-population allele/subtype frequencies.

### `imputationpanel/` — reference panel construction and evaluation

The largest pipeline: builds a phased Korea10K imputation reference panel and
benchmarks it against alternatives.

| Stage | Scripts |
| --- | --- |
| VCF/BCF preparation | `convert_pgen_to_bcf.py`, `convert_bcf_to_vcf.py`, `convert_chrname_vcf.py`, `make_biallelic_vcf_bcftools.py`, `make_vcf_index_bcftools.py`, `run_bgzip_vcf.py`, `run_tabix_vcf.py`, `fix_vcf_missing_contigs.py`, `multiple_run_fix_vcf.py` |
| Liftover | `run_liftover_gatk.py`, `clean_liftover_other_chromosome.py`, `run_selectvariants_gatk.py` |
| Chunking | `get_chr_chunks.py`, `add_header_chunk_file.py`, `split_chr_by_chunks_bcftools.py`, `split_chr_maf_info.py` |
| Phasing (SHAPEIT5) | `run_phase_common_shapeit.py`, `run_phase_common_shapeit_test.py`, `run_ligate_shapeit.py` |
| Panel build (Minimac) | `make_m3vcf_imputation_minimac3.py`, `make_msav_imputation_minimac4.py` |
| Imputation | `run_imputation_minimac.py` |
| Evaluation | `calc_concord_imp_glimpse.py`, `extract_R2_values_imputation_bcftools.py`, `prepare_input_correlation_imputation_genotype.py`, `calc_correlation_imputation_genotype.py`, `multiple_run_calc_corr_imp.py`, `concat_correlation_imputation_genotype.py`, `attach_maf_info_imputation_genotype.py`, `get_allele_frequency_plink.py`, `draw_imputation_performace_plot.py` |
| Chip sites | `parse_chip_manifest.py`, `get_samples_list_from_vcf.py` |

`hofmeister2023/` is **vendored third-party material** — the published scripts
and source data of Hofmeister et al. (2023), kept verbatim as the comparison
baseline. It is deliberately left unmodified, including its own paths.

### `methylation_bias/` — power loss from CpG-eliminating variants

A self-contained simulation study with its own
[README](analysis/methylation_bias/README.md), `requirements.txt` and
`run_all.sh`. It uses paths relative to its own directory and needs none of the
environment variables above.

### `variantdiscovery/` — variant discovery saturation

- `calc_variant_discovery_with_shuffling.py` / `plot_variant_discovery_with_shuffling.py`
  — discovery curves from repeatedly shuffled sample orderings, split by allele
  frequency class.
- `calc_allele_frequency.py` / `draw_allele_frequency.py` — allele frequency
  spectrum and the PlinkID↔rsID conversion tables it relies on.
- `Figure1.py` — a working draft of the corresponding manuscript panel.

---

## `paper/`

Scripts that render the manuscript figures and tables from the analysis outputs.
Each is standalone and writes into the figure/table output directory configured
through the environment variables.

| Script | Content |
| --- | --- |
| `Figure1.py` | Cohort overview — sample counts by category, major disease groups and disease subgroups |
| `Figure2.py` | Variant discovery — variant count against sample count and against alternate allele frequency |
| `Figure3.py` | Population structure — PCA of Korea10K with 1KGP; also exports supplementary Y-chromosome and mitochondrial tables |
| `Figure4.py` | CpG sites lost to sequence variation on methylation arrays, across 1KGP superpopulations and Korea10K |

| Script | Content |
| --- | --- |
| `Supplementary_Figures/Supplementary_Figure1.py` | Per-sample sequencing production and depth summary |
| `Supplementary_Figures/Supplementary_Figure2.py` | Sample counts by variant frequency category |
| `Supplementary_Figures/Supplementary_Figure3.py` | ADMIXTURE cross-validation error against *K* |
| `Supplementary_Figures/Supplementary_Figure7.py` | Shift of the methylation beta value away from the homozygous reference genotype at CpG-eliminating variants, by allele frequency |
| `Supplementary_Figures/Supplementary_Figure8.py` | Power of differential-methylation testing at CpG-eliminating variant sites (rendered from `analysis/methylation_bias/`) |
| `Supplementary_Figures/Supplementary_Figure9.py` | UpSet plot per chip — sharing of CpG-eliminating variants among AFR, AMR, EAS, EUR, SAS and Korea10K |
| `Supplementary_Figures/Supplementary_Figure10.py` | Distribution of average WGS depth across samples |
| `Supplementary_Tables/Supplementary_Table1.py` | Per-sample sequencing and depth metadata |
| `Supplementary_Tables/Supplementary_Table2.py` | Per-sample sequencing and depth metadata (companion breakdown) |
| `Supplementary_Tables/Supplementary_Table5.py` | Sequencing depth summarised by platform |

The supplementary numbering is not contiguous — the gaps correspond to items
that are not generated by code in this repository.

---

## Third-party tools

The scripts call out to the following, which are not bundled here:

PLINK 1.9 / 2.0 · ADMIXTURE 1.3 · SHAPEIT5 · Minimac3 / Minimac4 · GLIMPSE2 ·
bcftools · tabix / bgzip · GATK 4

Python dependencies are the usual scientific stack (`pandas`, `numpy`,
`matplotlib`, `seaborn`, `scipy`, `openpyxl`); `analysis/methylation_bias/`
pins its own set in `requirements.txt`.
