#!/usr/bin/env python3
"""
Publication figure: power loss from CpG-eliminating variants.

Emits PDF (vector, for submission), PNG (600 dpi), and the exact source data
behind every plotted point.

Two things differ from the exploratory version in 04_figures.py:

  * The carrier-free reference is drawn as the min-max band across the AF strata
    present in each panel, not as a single pooled curve. Pooling is misleading at
    n = 200, where carrier-free power itself ranges from 0.76 (rare variants) to
    0.42 (common ones) because the sites able to supply 200 intact samples are a
    biased subset at high AF. The band keeps every comparison AF-matched: a
    contaminated curve below the band has lost power against any stratum.

  * Monte Carlo error bars, so the reader can see the collapse is not noise.

--extended adds a second row with the genotype-adjusted model.
"""
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

mpl.use("Agg")
ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "figures"

AF_BINS = [0.0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 1.0]
AF_LAB = ["0-1%", "1-5%", "5-10%", "10-20%", "20-30%", "30-50%", ">50%"]
SEQ = LinearSegmentedColormap.from_list("seq", ["#9fc4ec", "#2a78d6", "#0d2647"])
AF_C = {l: SEQ(i / (len(AF_LAB) - 1)) for i, l in enumerate(AF_LAB)}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdcd8"
BAND = "#b8b8b2"
NS = [200, 400, 800]
ALPHA_BONF = 0.05 / 17338

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,          # embed as TrueType, editable
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.edgecolor": INK2, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "legend.frameon": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
})


def prepare():
    df = pd.read_parquet(RES / "sim_raw.parquet")
    df["af_bin"] = pd.cut(df.af, AF_BINS, labels=AF_LAB)
    real = df[(df.donor == "mix") & (df.frac < 0) & (df["mode"] == "random")]
    clean = df[df.donor == "none"]

    def agg(d, col):
        g = (d.groupby(["n", "af_bin", "delta"], observed=True)[col]
             .agg(["mean", "std", "size"]).reset_index())
        g["se"] = g["std"] / np.sqrt(g["size"])
        return g.rename(columns={"mean": col, "size": "n_sites"})

    return (agg(real, "rej_bonf"), agg(real, "rej_bonf_adj"),
            agg(clean, "rej_bonf").rename(columns={"rej_bonf": "clean"}))


def draw_row(axs, data, clean, value_col, show_xlabel):
    for ax, n in zip(axs, NS):
        cl = clean[clean.n == n]
        # AF-matched reference: the range of carrier-free power across the strata
        # that can actually reach this cohort size.
        band = cl.groupby("delta").clean.agg(["min", "max"]).sort_index()
        ax.fill_between(band.index, band["min"], band["max"], color=BAND,
                        alpha=0.45, lw=0, zorder=1)

        labels = []
        for lab in AF_LAB:
            g = data[(data.n == n) & (data.af_bin == lab)].sort_values("delta")
            if g.empty or g[value_col].isna().all():
                continue
            # 95% CI rather than +/-1 SE: at one SE the bars are the same size
            # as the markers and read as noise on the marker rather than as
            # uncertainty, which invites mistaking the grey band for the CI.
            ax.errorbar(g.delta, g[value_col], yerr=1.96 * g.se, color=AF_C[lab],
                        lw=1.6, marker="o", ms=3.0, capsize=2.4, elinewidth=1.1,
                        capthick=1.1, zorder=3)
            labels.append([g[value_col].iloc[-1], lab, AF_C[lab]])

        labels.sort(key=lambda r: r[0])
        for i in range(1, len(labels)):
            labels[i][0] = max(labels[i][0], labels[i - 1][0] + 0.062)
        for y, txt, c in labels:
            ax.text(0.207, y, txt, fontsize=6.8, color=c, va="center")

        ax.axhline(0.8, color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=2)
        ax.set_xlim(-0.008, 0.30)
        ax.set_ylim(-0.04, 1.30)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
        if show_xlabel:
            ax.set_xlabel("true methylation difference (Δβ)")
        ax.set_title(f"n = {n}", loc="left", fontweight="bold", pad=5)


def main(extended=False):
    real, real_adj, clean = prepare()
    nrow = 2 if extended else 1
    f, axs = plt.subplots(nrow, len(NS), figsize=(7.1, 2.5 * nrow + 0.5),
                          sharey=True, squeeze=False)

    draw_row(axs[0], real, clean, "rej_bonf", not extended)
    # name the test on the axis: the figure is read without the methods text
    axs[0][0].set_ylabel("power, Welch t-test\n(p < 2.9 × 10⁻⁶, Bonferroni)")
    if extended:
        draw_row(axs[1], real_adj, clean, "rej_bonf_adj", True)
        axs[1][0].set_ylabel("power, β ~ case + dosage\n(p < 2.9 × 10⁻⁶, Bonferroni)")
        axs[0][0].text(-0.30, 1.32, "A", transform=axs[0][0].transAxes,
                       fontsize=12, fontweight="bold")
        axs[1][0].text(-0.30, 1.32, "B", transform=axs[1][0].transAxes,
                       fontsize=12, fontweight="bold")

    axs[0][0].text(0.002, 0.835, "80% power", fontsize=6.5, color=INK2)
    handles = [
        Patch(facecolor=BAND, alpha=0.45,
              label="carrier-free reference: spread across AF strata "
                    "(not a confidence band)"),
        Line2D([], [], color=AF_C["1-5%"], marker="o", ms=3.0, lw=1.6,
               label="cohort with carriers, by allele frequency (error bars 95% CI)"),
    ]
    f.legend(handles=handles, fontsize=7, loc="lower center", ncol=2,
             bbox_to_anchor=(0.5, -0.035), handlelength=1.8, columnspacing=2.0)

    f.tight_layout(rect=(0, 0.015, 1, 0.99))
    for ext in ("pdf", "png"):
        out = FIG / (f"fig_power_extended.{ext}" if extended else f"fig_power.{ext}")
        f.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(f)

    # source data for every plotted point
    COLS = ["n", "af_bin", "delta", "power", "se", "n_sites"]

    def tidy(d, col, label):
        t = d[["n", "af_bin", "delta", col, "se", "n_sites"]].copy()
        t.columns = COLS
        t.insert(0, "series", label)
        return t

    out = pd.concat([
        tidy(clean, "clean", "no carriers (gt0-only cohort)"),
        tidy(real, "rej_bonf", "with carriers (HWE frequency, randomly split)"),
        tidy(real_adj, "rej_bonf_adj", "with carriers, genotype-adjusted model"),
    ])
    out.to_csv(RES / "fig_power_source_data.tsv", sep="\t", index=False,
               float_format="%.4f")

    tag = "fig_power_extended" if extended else "fig_power"
    print(f"wrote figures/{tag}.pdf, figures/{tag}.png, "
          f"results/fig_power_source_data.tsv")
    print(f"max Monte Carlo SE across plotted cells: {real.se.max():.4f}"
          f"  (95% CI half-width {1.96 * real.se.max():.4f})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--extended", action="store_true",
                    help="add a second row with the genotype-adjusted model")
    a = ap.parse_args()
    main(a.extended)
