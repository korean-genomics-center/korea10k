#!/usr/bin/env python3
"""
Spike-in permutation simulation: does a CpG-eliminating variant bias DMS calling?

DESIGN
------
Synthetic case/control cohorts are built out of each site's own observed beta
values, so the noise structure is real rather than parametric.

  cohort      n in {200,400,600,800,1000} drawn without replacement from the
              site's gt0 (intact-CpG) samples, split 5:5 into case / control.
  spike-in    a shift delta in {0,.01,.05,.10,.15,.20} added to the case group,
              signed toward 0.5 (down at hypermethylated sites, up at
              hypomethylated ones) so it never runs into the 0/1 boundary.
  donors      a fraction f of the cohort is replaced by variant carriers:
                gt2  hom-alt, CpG destroyed on both alleles
                gt1  het, the far more common carrier class at most AF
                mix  gt1 and gt2 in the HWE ratio 2p(1-p) : p^2
  assignment  how the donors distribute across the two arms:
                balanced  exactly ktot/2 in each arm  (perfect genotype matching)
                random    Binomial(ktot, 0.5) into the case arm -- what a real
                          cohort looks like when genotype is independent of
                          phenotype, because carrier counts differ by chance
                caseonly  all donors in the case arm -- genotype fully
                          confounded with case status
  tests       (1) naive Welch t-test on beta, case vs control
              (2) OLS  beta ~ case + dosage   (genotype-adjusted; the mitigation)

WHY 'random' MATTERS: 'balanced' fixes the carrier count per arm, which removes
the between-genotype component from the true sampling variance of the mean
difference while leaving it in the variance the t-test estimates. That makes the
test anti-conservative-looking in reverse -- it is too conservative, an artefact
of the matching rather than a property of the data. 'random' is the honest null.

MODELLING CHOICE THAT MATTERS
-----------------------------
The spike-in is applied only to samples whose CpG is intact. A variant that
destroys the CpG removes the substrate, so a disease-associated methylation
change cannot occur there. Contamination therefore both inflates variance and
*dilutes* the true effect. Use --spike-all to shift every sample instead.

Because the shift is additive and boundary-free, every delta is obtained in
closed form from the per-replicate sufficient statistics of a single draw, which
is what makes the full grid tractable.

Output: results/sim_raw.parquet, one row per
        (site, n, donor, mode, frac, delta), aggregated over --reps replicates.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent

NS = np.array([200, 400, 600, 800, 1000])
DELTAS = np.array([0.0, 0.01, 0.05, 0.10, 0.15, 0.20])
# HWE_F is a sentinel meaning "use this site's own carrier frequency"
HWE_F = -1.0
FRACS = np.array([0.01, 0.025, 0.05, 0.10, 0.20, 0.35, HWE_F])
MODES = ("balanced", "random", "caseonly")
DONORS = ("gt2", "gt1", "mix")
ALPHA_NOM = 0.05
N_SITES_TOTAL = 17338
ALPHA_BONF = 0.05 / N_SITES_TOTAL


def welch(sum_a, sq_a, na, sum_b, sq_b, nb):
    """Welch t-test from sufficient statistics; arrays broadcast elementwise."""
    ma, mb = sum_a / na, sum_b / nb
    va = np.maximum((sq_a - na * ma * ma) / (na - 1), 1e-12)
    vb = np.maximum((sq_b - nb * mb * mb) / (nb - 1), 1e-12)
    sa, sb = va / na, vb / nb
    se = np.sqrt(sa + sb)
    df = (sa + sb) ** 2 / (sa * sa / (na - 1) + sb * sb / (nb - 1))
    return ma - mb, se, (ma - mb) / se, df, va, vb


def ols_adjusted(N, s_case, s_dos, s_casedos, s_dos2, s_y, s_casey, s_dosy, s_y2):
    """
    OLS of y on [1, case, dosage]; returns (beta_case, se, t, df, collinear).
    Solved from cross-product sufficient statistics so the whole
    replicate x delta grid is done in one vectorised pass.
    """
    mc, md, my = s_case / N, s_dos / N, s_y / N
    scc = s_case - N * mc * mc          # case is a 0/1 indicator
    sdd = s_dos2 - N * md * md
    scd = s_casedos - N * mc * md
    scy = s_casey - N * mc * my
    sdy = s_dosy - N * md * my
    syy = s_y2 - N * my * my

    det = scc * sdd - scd * scd
    # tol scaled to the design: det vanishes when dosage is constant (no donors)
    # or perfectly collinear with case status (every case sample is a carrier).
    tol = 1e-9 * np.maximum(scc * sdd, 1e-12)
    degen = det <= tol
    det_s = np.where(degen, 1.0, det)
    scc_s = np.maximum(scc, 1e-12)

    b_case = np.where(degen, scy / scc_s, (sdd * scy - scd * sdy) / det_s)
    b_dos = np.where(degen, 0.0, (scc * sdy - scd * scy) / det_s)

    rss = np.maximum(syy - b_case * scy - b_dos * sdy, 1e-12)
    df = np.where(degen, N - 2, N - 3)
    sigma2 = rss / df
    var_b = np.where(degen, sigma2 / scc_s, sigma2 * sdd / det_s)
    se = np.sqrt(np.maximum(var_b, 1e-24))
    return b_case, se, b_case / se, df, degen


def cumstats(*arrays):
    """Cumulative sums along axis 1, each with a leading zero column."""
    z = np.zeros((arrays[0].shape[0], 1))
    return [np.concatenate([z, np.cumsum(a, axis=1)], axis=1) for a in arrays]


def take(cum, k):
    """cum[:, k] where k is a scalar or a per-row index array."""
    if np.isscalar(k) or np.ndim(k) == 0:
        return cum[:, int(k)]
    return cum[np.arange(cum.shape[0]), k]


def run(args):
    d = np.load(ROOT / "data" / "betas.npz", allow_pickle=False)
    off0, off1, off2 = d["off0"], d["off1"], d["off2"]
    v0 = d["v0"].astype(np.float64)
    v1 = d["v1"].astype(np.float64)
    v2 = d["v2"].astype(np.float64)
    af, chrom, pos = d["af"], d["chrom"], d["pos"]

    sites = np.load(ROOT / "data" / "site_selection.npy")
    rng = np.random.default_rng(args.seed)
    R = args.reps
    zeroR = np.zeros(R)

    rows = []
    t0 = time.time()
    for si, s in enumerate(sites):
        b0 = v0[off0[s]:off0[s + 1]]
        g1 = v1[off1[s]:off1[s + 1]]
        g2 = v2[off2[s]:off2[s + 1]]
        p = float(af[s])
        f_hwe = 2 * p * (1 - p) + p * p            # carrier frequency under HWE
        sgn = -1.0 if b0.mean() > 0.5 else 1.0     # shift toward 0.5, stays on-scale
        a = sgn * DELTAS
        meta = dict(site=int(s), chrom=str(chrom[s]), pos=int(pos[s]),
                    af=p, sgn=sgn, mean_gt0=float(b0.mean()), f_hwe=f_hwe)

        # donor pools as (value, dosage) pairs
        pools = {}
        if len(g2):
            pools["gt2"] = (g2, np.full(len(g2), 2.0))
        if len(g1):
            pools["gt1"] = (g1, np.full(len(g1), 1.0))
        if len(g1) and len(g2):
            pools["mix"] = (np.concatenate([g1, g2]),
                            np.concatenate([np.ones(len(g1)), np.full(len(g2), 2.0)]))
        # HWE sampling weights for the mix pool: het and hom classes enter at
        # 2p(1-p) : p^2, uniform within each class.
        mix_w = None
        if "mix" in pools:
            w = np.concatenate([np.full(len(g1), 2 * p * (1 - p) / len(g1)),
                                np.full(len(g2), p * p / len(g2))])
            mix_w = w / w.sum()

        ns_here = NS[NS <= len(b0)]
        if len(ns_here) == 0:
            continue

        # one permutation of the gt0 pool per replicate; the first n entries are
        # a uniform without-replacement sample, so cohort sizes are nested draws.
        perm = np.argsort(rng.random((R, len(b0))), axis=1)

        for n in ns_here:
            h = n // 2
            idx = perm[:, :n]
            (cs_c,), (cq_c,) = cumstats(b0[idx[:, :h]]), cumstats(b0[idx[:, :h]] ** 2)
            (cs_k,), (cq_k,) = cumstats(b0[idx[:, h:]]), cumstats(b0[idx[:, h:]] ** 2)
            T_s_c, T_q_c = cs_c[:, -1], cq_c[:, -1]
            T_s_k, T_q_k = cs_k[:, -1], cq_k[:, -1]

            # ---------- clean arm (no carriers at all) ----------
            emit(rows, meta, n, "none", "none", 0.0, 0, 0, False, a, args.spike_all,
                 T_s_c, T_q_c, T_s_k, T_q_k,
                 zeroR, zeroR, zeroR, zeroR, zeroR, h, h,
                 intact_s_c=T_s_c, n_intact_c=h)

            # ---------- contaminated arms ----------
            for donor, (dvals_pool, ddose_pool) in pools.items():
                fr = np.array([f if f != HWE_F else f_hwe for f in FRACS])
                kmax = int(min(h * 2, round(fr.max() * n)))
                if kmax < 1:
                    continue
                wr = len(dvals_pool) < kmax
                if wr:
                    sel = (rng.choice(len(dvals_pool), size=(R, kmax), p=mix_w)
                           if donor == "mix"
                           else rng.integers(0, len(dvals_pool), size=(R, kmax)))
                else:
                    key = rng.random((R, len(dvals_pool)))
                    if donor == "mix":                       # weighted w/o replacement
                        key = np.log(np.maximum(key, 1e-300)) / mix_w[None, :]
                        sel = np.argsort(-key, axis=1)[:, :kmax]
                    else:
                        sel = np.argsort(key, axis=1)[:, :kmax]
                dv_, dd_ = dvals_pool[sel], ddose_pool[sel]
                ds, dq, dD, dD2, dDY = cumstats(dv_, dv_ ** 2, dd_, dd_ ** 2, dd_ * dv_)

                for f, f_nom in zip(fr, FRACS):
                    ktot = int(round(f * n))
                    if ktot < 1:
                        continue
                    for mode in MODES:
                        if mode == "balanced":
                            kc = np.full(R, ktot // 2)
                        elif mode == "random":
                            kc = rng.binomial(ktot, 0.5, size=R)
                        else:
                            kc = np.full(R, ktot)
                        kk = ktot - kc
                        if kc.max() > h or kk.max() > h:
                            continue
                        # donors occupy columns [0,kc) of case and [kc,kc+kk) of control
                        d_s_c, d_q_c = take(ds, kc), take(dq, kc)
                        d_D_c, d_D2_c, d_DY_c = take(dD, kc), take(dD2, kc), take(dDY, kc)
                        e = kc + kk
                        d_s_k, d_q_k = take(ds, e) - d_s_c, take(dq, e) - d_q_c
                        d_D_k = take(dD, e) - d_D_c
                        d_D2_k = take(dD2, e) - d_D2_c
                        d_DY_k = take(dDY, e) - d_DY_c
                        # intact-CpG members surviving the substitution
                        g_s_c, g_q_c = T_s_c - take(cs_c, kc), T_q_c - take(cq_c, kc)
                        g_s_k, g_q_k = T_s_k - take(cs_k, kk), T_q_k - take(cq_k, kk)
                        emit(rows, meta, n, donor, mode, f_nom, kc, kk, wr, a,
                             args.spike_all,
                             g_s_c + d_s_c, g_q_c + d_q_c, g_s_k + d_s_k, g_q_k + d_q_k,
                             d_D_c, d_D_k, d_D2_c + d_D2_k, d_DY_c, d_DY_k,
                             h, h, intact_s_c=g_s_c, n_intact_c=h - kc)

        if args.progress and (si + 1) % 100 == 0:
            el = time.time() - t0
            print(f"  {si+1}/{len(sites)} sites  {el:.0f}s  "
                  f"eta {el/(si+1)*(len(sites)-si-1):.0f}s", file=sys.stderr, flush=True)

    df = pd.DataFrame(rows)
    out = ROOT / "results" / args.out
    df.to_parquet(out, index=False)
    print(f"wrote {out}  rows={len(df)}  ({time.time()-t0:.0f}s)", file=sys.stderr)


def emit(rows, meta, n, donor, mode, frac, kc, kk, wr, a, spike_all,
         S_c, Q_c, S_k, Q_k, D_c, D_k, D2, DY_c, DY_k,
         nc, nk, intact_s_c, n_intact_c):
    """
    Apply every spike-in level in closed form, run both tests, append one
    summary row per delta. S_*/Q_* are the (R,) pre-spike sum / sum-of-squares.
    """
    col = lambda v: np.asarray(v, dtype=float).reshape(-1, 1)
    m = col(nc if spike_all else n_intact_c)        # case samples that shift
    src = col(S_c if spike_all else intact_s_c)     # their pre-spike sum

    Sc = col(S_c) + a[None, :] * m
    Qc = col(Q_c) + 2 * a[None, :] * src + (a * a)[None, :] * m
    ones = np.ones_like(Sc)
    Sk, Qk = col(S_k) * ones, col(Q_k) * ones

    diff, se, t, df_, va, vb = welch(Sc, Qc, nc, Sk, Qk, nk)
    p = 2 * stats.t.sf(np.abs(t), df_)

    N = nc + nk
    # dosage . y : carriers never receive the spike unless --spike-all
    dosy = col(DY_c + DY_k) * ones
    if spike_all:
        dosy = dosy + a[None, :] * col(D_c)
    b_adj, se_adj, t_adj, dfa, degen = ols_adjusted(
        N, float(nc), col(D_c + D_k) * ones, col(D_c) * ones, col(D2) * ones,
        Sc + Sk, Sc, dosy, Qc + Qk)
    p_adj = 2 * stats.t.sf(np.abs(t_adj), dfa)

    sgn = meta["sgn"]
    sdiff, sadj = sgn * diff, sgn * b_adj           # sign so a clean run gives +delta
    base = {k: v for k, v in meta.items() if k != "sgn"}
    for j, delta in enumerate(DELTAS):
        rows.append(dict(
            **base, n=int(n), donor=donor, mode=mode, frac=float(frac),
            k_case=float(np.mean(kc)), k_ctrl=float(np.mean(kk)),
            delta=float(delta), with_replacement=bool(wr),
            collinear=float(np.mean(degen[:, j])),
            dhat=float(sdiff[:, j].mean()), dhat_sd=float(sdiff[:, j].std(ddof=1)),
            se_mean=float(se[:, j].mean()),
            var_case=float(va[:, j].mean()), var_ctrl=float(vb[:, j].mean()),
            rej05=float((p[:, j] < ALPHA_NOM).mean()),
            rej_bonf=float((p[:, j] < ALPHA_BONF).mean()),
            dhat_adj=float(sadj[:, j].mean()), dhat_adj_sd=float(sadj[:, j].std(ddof=1)),
            se_adj_mean=float(se_adj[:, j].mean()),
            rej05_adj=float((p_adj[:, j] < ALPHA_NOM).mean()),
            rej_bonf_adj=float((p_adj[:, j] < ALPHA_BONF).mean()),
        ))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20240301)
    ap.add_argument("--spike-all", action="store_true")
    ap.add_argument("--out", default="sim_raw.parquet")
    ap.add_argument("--progress", action="store_true")
    run(ap.parse_args())
