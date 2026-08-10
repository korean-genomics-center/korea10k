#!/usr/bin/env python3
"""
Check the closed-form machinery in 02_simulate.py against a brute-force
simulation that materialises every cohort explicitly.

Verified on identical random draws:
  1. the sufficient-statistic expansion of the spike-in reproduces the Welch
     test computed on explicit spiked arrays;
  2. ols_adjusted reproduces statsmodels OLS of beta ~ case + dosage, including
     the mixed-dosage (gt1 + gt2) case and per-replicate varying carrier counts.
"""
import sys
from importlib import import_module
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
sim = import_module("02_simulate")

ROOT = Path(__file__).resolve().parent.parent


def draws(b0, dpool, ddose, n, ktot, mode, reps, seed):
    """Generate the shared cohort draws both paths consume."""
    g = np.random.default_rng(seed)
    h = n // 2
    out = []
    for _ in range(reps):
        idx = g.permutation(len(b0))[:n]
        if mode == "balanced":
            kc = ktot // 2
        elif mode == "random":
            kc = int(g.binomial(ktot, 0.5))
        else:
            kc = ktot
        kk = ktot - kc
        sel = g.choice(len(dpool), ktot, replace=len(dpool) < ktot)
        out.append((idx, kc, kk, sel))
    return h, out


def brute(b0, dpool, ddose, n, ktot, mode, delta, sgn, reps, seed):
    h, dr = draws(b0, dpool, ddose, n, ktot, mode, reps, seed)
    a = sgn * delta
    diffs, ps, badj, padj = [], [], [], []
    for idx, kc, kk, sel in dr:
        case, ctrl = b0[idx[:h]].copy(), b0[idx[h:]].copy()
        dose_c, dose_k = np.zeros(h), np.zeros(h)
        if kc:
            case[:kc] = dpool[sel[:kc]]
            dose_c[:kc] = ddose[sel[:kc]]
        if kk:
            ctrl[:kk] = dpool[sel[kc:kc + kk]]
            dose_k[:kk] = ddose[sel[kc:kc + kk]]
        case[dose_c == 0] += a                      # spike intact CpGs only
        tt = stats.ttest_ind(case, ctrl, equal_var=False)
        diffs.append(case.mean() - ctrl.mean())
        ps.append(tt.pvalue)
        y = np.concatenate([case, ctrl])
        cs = np.concatenate([np.ones(h), np.zeros(h)])
        dose = np.concatenate([dose_c, dose_k])
        X = sm.add_constant(cs if np.ptp(dose) == 0
                            else np.column_stack([cs, dose]))
        fit = sm.OLS(y, X).fit()
        badj.append(fit.params[1])
        padj.append(fit.pvalues[1])
    return (sgn * np.mean(diffs), np.mean(np.array(ps) < 0.05),
            sgn * np.mean(badj), np.mean(np.array(padj) < 0.05))


def fast(b0, dpool, ddose, n, ktot, mode, delta, sgn, reps, seed):
    h, dr = draws(b0, dpool, ddose, n, ktot, mode, reps, seed)
    R = reps
    Sc = np.empty(R); Qc = np.empty(R); Sk = np.empty(R); Qk = np.empty(R)
    Dc = np.empty(R); Dk = np.empty(R); D2 = np.empty(R)
    DYc = np.empty(R); DYk = np.empty(R)
    dsc = np.empty(R); dsk = np.empty(R); isc = np.empty(R); nic = np.empty(R)
    for i, (idx, kc, kk, sel) in enumerate(dr):
        case, ctrl = b0[idx[:h]].copy(), b0[idx[h:]].copy()
        dose_c, dose_k = np.zeros(h), np.zeros(h)
        if kc:
            case[:kc] = dpool[sel[:kc]]
            dose_c[:kc] = ddose[sel[:kc]]
        if kk:
            ctrl[:kk] = dpool[sel[kc:kc + kk]]
            dose_k[:kk] = ddose[sel[kc:kc + kk]]
        Sc[i], Qc[i] = case.sum(), (case ** 2).sum()
        Sk[i], Qk[i] = ctrl.sum(), (ctrl ** 2).sum()
        Dc[i], Dk[i] = dose_c.sum(), dose_k.sum()
        D2[i] = (dose_c ** 2).sum() + (dose_k ** 2).sum()
        DYc[i], DYk[i] = (dose_c * case).sum(), (dose_k * ctrl).sum()
        dsc[i], dsk[i] = case[:kc].sum(), ctrl[:kk].sum()
        isc[i], nic[i] = case[kc:].sum(), h - kc
    rows = []
    sim.DELTAS = np.array([delta])
    meta = dict(site=0, chrom="chr1", pos=1, af=0.1, sgn=sgn,
                mean_gt0=float(b0.mean()), f_hwe=0.1)
    sim.emit(rows, meta, n, "d", mode, 0.1, dsc * 0, dsk * 0, False,
             np.array([sgn * delta]), False,
             Sc, Qc, Sk, Qk, Dc, Dk, D2, DYc, DYk, h, h,
             intact_s_c=isc, n_intact_c=nic)
    r = rows[0]
    return r["dhat"], r["rej05"], r["dhat_adj"], r["rej05_adj"]


def main():
    d = np.load(ROOT / "data" / "betas.npz", allow_pickle=False)
    off0, off1, off2 = d["off0"], d["off1"], d["off2"]
    v0, v1, v2 = (d[k].astype(float) for k in ("v0", "v1", "v2"))
    sites = np.load(ROOT / "data" / "site_selection.npy")
    cand = [s for s in sites
            if d["n0"][s] >= 400 and d["n2"][s] >= 60 and d["n1"][s] >= 60][:3]

    print(f"{'site':>6} {'donor':>5} {'mode':>9} {'delta':>6} | "
          f"{'dhat brute':>10} {'dhat fast':>10} | {'rej b':>6} {'rej f':>6} | "
          f"{'adj brute':>10} {'adj fast':>10} | {'arej b':>6} {'arej f':>6}")
    worst = 0.0
    for s in cand:
        b0 = v0[off0[s]:off0[s + 1]]
        g1 = v1[off1[s]:off1[s + 1]]
        g2 = v2[off2[s]:off2[s + 1]]
        sgn = -1.0 if b0.mean() > 0.5 else 1.0
        pools = {
            "gt2": (g2, np.full(len(g2), 2.0)),
            "mix": (np.concatenate([g1, g2]),
                    np.concatenate([np.ones(len(g1)), np.full(len(g2), 2.0)])),
        }
        for dn, (dp, dd) in pools.items():
            for mode in ("balanced", "random", "caseonly"):
                for delta in (0.0, 0.10):
                    args = (b0, dp, dd, 400, 40, mode, delta, sgn, 150, 4242)
                    bb, ff = brute(*args), fast(*args)
                    worst = max(worst, max(abs(x - y) for x, y in zip(bb, ff)))
                    print(f"{s:>6} {dn:>5} {mode:>9} {delta:>6.2f} | "
                          f"{bb[0]:>10.5f} {ff[0]:>10.5f} | {bb[1]:>6.2f} {ff[1]:>6.2f} | "
                          f"{bb[2]:>10.5f} {ff[2]:>10.5f} | {bb[3]:>6.2f} {ff[3]:>6.2f}")
    print(f"\nmax absolute discrepancy: {worst:.3e}")
    print("PASS" if worst < 1e-9 else "FAIL")
    sys.exit(0 if worst < 1e-9 else 1)


if __name__ == "__main__":
    main()
