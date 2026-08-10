# Power loss from CpG-eliminating variants in differential methylation testing

Code to reproduce **Supplementary Figure S*x*** of *[paper title]* ([DOI]).

A single-nucleotide variant that disrupts a CpG dinucleotide removes the
substrate for methylation on the allele carrying it, so the beta value at such a
site reports genotype as much as biology. This repository contains the
simulation that quantifies what that does to the power of a conventional
differential-methylation analysis, and whether adjusting for local genotype
recovers it.

---

## What the figure shows

Case–control cohorts are assembled from the *observed* beta values at each
CpG-eliminating variant site, a true methylation difference Δβ is spiked into the
case arm, and the proportion of replicates reaching a Bonferroni-corrected
threshold is recorded.

- **Panel A** — conventional analysis (Welch two-sample *t*-test on beta).
- **Panel B** — the same cohorts with genotype adjustment (`beta ~ case + dosage`).
- **Coloured lines** — cohorts containing variant carriers at each site's own
  Hardy–Weinberg frequency, allocated between arms at random, so genotype is
  *independent of case status*. Shaded by allele frequency.
- **Grey band** — the carrier-free reference, drawn as the range across allele
  frequency strata so that comparisons are frequency-matched. It is a spread
  across strata, **not** a confidence interval.

Because the simulation uses observed values rather than a parametric model, the
skewness of the beta distribution, the pile-up at the 0 and 1 boundaries, and
the genotype-specific heteroscedasticity are all preserved.

---

## Repository contents

```
scripts/
  01_parse.py                 input TSV -> compact binary store (data/betas.npz)
  00_select_sites.py          stratified site selection by allele frequency
  02_simulate.py              the simulation engine
  02b_validate.py             checks the fast path against a brute-force implementation
  06_figure_power.py          builds the figure and its source-data table
  make_synthetic_input.py     generates a format-compatible synthetic dataset
```

Script numbering follows the full analysis for the paper; only the subset needed
for this figure is released here, which is why the numbers are not contiguous.
The five analysis scripts are byte-identical to those used to produce the
published figure, with one exception noted under *Adapting to other data*.

---

## Requirements

Python 3.11+ and the packages in `requirements.txt`. Versions used for the
published run:

```
numpy 2.3.5   pandas 2.3.3   scipy 1.16.3
statsmodels 0.14.5   matplotlib 3.10.6   pyarrow 21.0.0
```

```bash
pip install -r requirements.txt
```

`statsmodels` is needed only by the validation script.

---

## Running it

### With the study data

Place the input file (see `data/README.md` for the schema and how to obtain it)
in the repository root as
`per_cpg_beta_by_genotype.with_values.filtered.tsv`, then:

```bash
python3 scripts/01_parse.py            # ~1 min
python3 scripts/00_select_sites.py     # seconds
python3 scripts/02b_validate.py        # ~1 min, must print PASS
python3 scripts/02_simulate.py --progress   # ~2-3 min, writes ~120 MB
python3 scripts/06_figure_power.py --extended
```

or simply `bash run_all.sh`. Output:

```
figures/fig_power_extended.pdf      the figure, vector
figures/fig_power_extended.png      the figure, 600 dpi
results/fig_power_source_data.tsv   the value behind every plotted point
results/sim_raw.parquet             the full simulation grid
```

Omit `--extended` for a single-row version with panel A only.

### Without the study data

The input contains individual-level methylation values and is not distributed
here. To exercise the code path, generate a synthetic dataset with the same
schema:

```bash
python3 scripts/make_synthetic_input.py --sites 2000
bash run_all.sh
```

This reproduces the *shape* of the analysis but not its results: the values come
from a parametric model and carry no biological information. The figure it
produces will not match the published one, and should not be interpreted.

---

## How the simulation works

For each site, each of *R* = 100 replicates:

1. Draw *n* samples without replacement from the site's intact-CpG (gt0) pool
   and split them equally into case and control arms. Cohort sizes 200, 400,
   600, 800 and 1,000 are simulated; the figure shows 200, 400 and 800.
2. Substitute variant carriers in at the site's Hardy–Weinberg frequency
   *f* = 1 − (1 − *p*)², drawn from that site's own observed heterozygotes and
   alternate homozygotes in the ratio 2*p*(1 − *p*) : *p*². The carrier fraction
   is fixed by allele frequency, not chosen.
3. Allocate carriers between arms as Binomial(*k*, 0.5), so genotype is
   independent of case status while realised counts differ by chance, as in a
   real cohort.
4. Add Δβ to the case arm, signed toward 0.5 so it never reaches the 0 or 1
   boundary, and applied only to samples whose CpG is intact — a destroyed CpG
   cannot carry a methylation change. Contamination therefore both inflates
   variance and dilutes the true effect.
5. Test with Welch's *t*-test and with `beta ~ case + dosage`.

Power is the proportion of replicates rejected at α = 0.05/17,338 = 2.88 × 10⁻⁶,
averaged over sites within an allele-frequency stratum.

**On the exactly-balanced alternative.** Substituting an equal number of
carriers into each arm is *not* a valid null. Matching the arms on genotype
removes the between-genotype component from the true sampling variance of the
mean difference while leaving it in the variance the *t*-test estimates, so the
test rejects at only 0.003–0.026 against a nominal 0.05 and appears harmless.
`02_simulate.py` evaluates that allocation alongside the others — its rows carry
`mode == "balanced"` in `results/sim_raw.parquet`, so the comparison can be
reproduced from the same output — while the figure uses the binomial allocation.

**Why it is fast.** The spike-in is additive and never truncated, so every Δβ is
obtained in closed form from the sufficient statistics of a single draw rather
than by re-simulating — exactly, not approximately. `02b_validate.py` verifies
this against a brute-force implementation that materialises every cohort and
calls `scipy.stats.ttest_ind` and `statsmodels` OLS; the maximum discrepancy is
6 × 10⁻¹⁶, at floating-point precision. **Run it before trusting any output.**
As a second check, carrier-free cohorts reject at 0.043–0.053 across all
conditions, confirming the machinery is calibrated.

`02_simulate.py` computes a superset of what the figure uses — additional
contamination levels, single-genotype carrier classes, and allocations in which
genotype is correlated with case status. The figure selects
`donor == "mix"`, `frac < 0` (the Hardy–Weinberg frequency) and
`mode == "random"`, plus the carrier-free arm. The other conditions support
results reported elsewhere in the paper and are left in place so that the
released code is exactly what was run.

---

## Determinism

Random number generation uses NumPy's PCG64 with fixed seeds (20240301 for the
simulation, 11 for site selection, 4242 for validation). Re-running on the same
input reproduces the figure bit for bit. Results are not guaranteed identical
across NumPy major versions, since the generator's stream is version-stable but
the reduction order in aggregation is not.

---

## Adapting to other data

Two constants are specific to this study and must be changed if you apply this
to a different panel:

- `N_SITES_TOTAL = 17338` in `02_simulate.py` is the multiple-testing
  denominator, the size of the full variant panel — *not* the number of
  simulated sites.
- `AF_BINS` / `PER_BIN` in `00_select_sites.py` set the stratification.

`01_parse.py` takes an optional `--input` argument in this release, so the file
need not be renamed; the internal version had the path hard-coded. This is the
only difference from the code that produced the published figure.

---

## Runtime

On a 2023 Apple M2 (10 cores), with 1,540 sites and 100 replicates: parsing
~60 s, simulation ~135 s, figure ~10 s. Peak memory ~4 GB. The simulation writes
a ~120 MB parquet file.

---

## Data availability

The input is individual-level DNA methylation and genotype data and is **not**
included in this repository. See `data/README.md` for the file schema and for
access conditions.

---

## Citation

```bibtex
@article{[key],
  title   = {[title]},
  author  = {[authors]},
  journal = {[journal]},
  year    = {[year]},
  doi     = {[doi]}
}
```

## License

[Choose one — MIT is a common default for analysis code; see LICENSE.]
