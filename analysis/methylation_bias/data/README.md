# Input data

## Availability

The analysis takes per-individual DNA methylation beta values, grouped by
genotype, at CpG sites disrupted by a nearby single-nucleotide variant. **These
are individual-level human data and are not distributed in this repository.**

Access conditions: *[data access statement — repository/accession, controlled
access committee, or contact]*

`scripts/make_synthetic_input.py` generates a file with the same schema so that
the code can be run without these data. It is a format demonstration only; the
values come from a parametric model and results obtained from it will not match
the published figure.

> **If you fork this repository, do not commit the input file.** It is listed in
> `.gitignore` for that reason. Per-individual beta values at variant sites can
> carry genotype information, and a public copy would be a disclosure of
> individual-level data regardless of how the identifiers are labelled.

## Schema

A tab-separated file, one row per CpG site, with a header line and 17 columns.

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `chr` | string | chromosome, e.g. `chr1` |
| 2 | `pos` | integer | 1-based position of the CpG |
| 3 | `n_gt0` | integer | number of individuals with genotype 0/0 (CpG intact) |
| 4 | `mean_gt0` | float | mean beta among those individuals |
| 5 | `var_gt0` | float | variance of beta among those individuals |
| 6 | `med_gt0` | float | median beta among those individuals |
| 7–10 | `n_gt1`, `mean_gt1`, `var_gt1`, `med_gt1` | | as above, genotype 0/1 |
| 11–14 | `n_gt2`, `mean_gt2`, `var_gt2`, `med_gt2` | | as above, genotype 1/1 |
| 15 | `beta_gt0` | string | comma-separated beta values, one per 0/0 individual |
| 16 | `beta_gt1` | string | comma-separated beta values, one per 0/1 individual |
| 17 | `beta_gt2` | string | comma-separated beta values, one per 1/1 individual |

Requirements the pipeline relies on:

- Columns 15–17 must contain exactly `n_gt0`, `n_gt1` and `n_gt2` values
  respectively. `01_parse.py` will fail loudly on a mismatch.
- Beta values are in [0, 1]. The study file carries four decimal places.
- Individuals need not be in a consistent order across rows; the simulation
  treats each site independently and never joins values across sites.
- Genotype 1 is the allele that eliminates the CpG, so beta is expected to fall
  as genotype goes 0/0 → 0/1 → 1/1. Nothing enforces this, but the analysis is
  only meaningful for sites where it holds.
- Sites with an empty genotype class are permitted in the file;
  `00_select_sites.py` excludes them, since the simulation needs at least one
  carrier of each alternate genotype.

Allele frequency is derived, not read from the file:
`p = (n_gt1 + 2 * n_gt2) / (2 * (n_gt0 + n_gt1 + n_gt2))`.

### Example

```
chr	pos	n_gt0	mean_gt0	var_gt0	med_gt0	n_gt1	...	beta_gt0	beta_gt1	beta_gt2
chr1	903351	1120	0.95566	0.00219	0.96296	82	...	1.0000,0.9130,...	0.5806,...	0.0000,0.0000
```

## Study dataset

The published figure used 17,338 CpG-eliminating variants with 1,201–1,213
individuals per site (median 1,209), from the Korea10K cohort. The
multiple-testing denominator in `02_simulate.py` (`N_SITES_TOTAL = 17338`) refers
to this panel size and should be changed for any other dataset.

## Generated files

`01_parse.py` and `00_select_sites.py` write two intermediates here. Both are
derived and are excluded from version control.

| File | Written by | Contents |
|---|---|---|
| `betas.npz` | `01_parse.py` | beta values as flat float32 arrays with per-site offsets, plus genotype counts and allele frequency |
| `site_selection.npy` | `00_select_sites.py` | integer indices of the sites the simulation runs on |
