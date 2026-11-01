# -*- coding: utf-8 -*-
"""run_extended.py -- extended analyses for the Dual-Lock manuscript (v3 of the study).

    python oral_peptide_mc.py --test ; python oral_peptide_mc.py --run     # original model and v1 results (unchanged)
    python run_extended.py --test                                           # verification tests of the extension
    python run_extended.py [stage ...]                                      # stages -> results_ext/ext_<stage>.json (default: all)

Stages: equiv paired converge corr sobol pareto fals struct abc valid
All draws use fixed seeds; every arm of a comparison is run on the same uniforms (common random numbers).
"""
import json
import os
import sys
import time

import numpy as np
from scipy.stats import rankdata, spearmanr

import oral_peptide_mc as M
import dualock_lib as L

OUT = "results_ext"
S0 = 20261101
N_MAIN = 20000
P0_LUMPED = None          # calibrated at start (same procedure and value as v1: 0.009746...)
ARM_DEF = M.ARMS
CORR_BASE = {("EC50", "C_tox"): 0.6, ("EC50", "Emax"): -0.3, ("C_tox", "t_rev"): -0.4}   # grade-D assumptions (see manuscript Table 2)
T0 = time.time()


def tick(msg):
    print("[%6.0fs] %s" % (time.time() - T0, msg), flush=True)


def save(name, obj):
    os.makedirs(OUT, exist_ok=True)
    json.dump(obj, open(os.path.join(OUT, "ext_%s.json" % name), "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


def p0():
    global P0_LUMPED
    if P0_LUMPED is None:
        P0_LUMPED = M.calibrate()
    return P0_LUMPED


def run(spec, U, **kw):
    return L.run_u(U, spec, p0(), **kw)


def paired_ci(a, b, z=1.96):
    d = a.astype(float) - b.astype(float)
    m = d.mean(); se = d.std(ddof=1) / np.sqrt(len(d))
    return [float(m), float(m - z * se), float(m + z * se)]


# ------------------------------------------------------------------------------------------------ stages
def stage_equiv():
    """The U-based sampler is a re-implementation of the v1 sampler: same marginals, different random-number order. Check the headline arms."""
    v1 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "results_v1.json"))) if os.path.exists(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "results_v1.json")) else None
    U = L.uniforms(N_MAIN, S0 + 1)
    out = {}
    for arm in ("R1_reference_PE", "P1_dual_lock", "P3_no_colocalisation", "P5_no_self_limiting_PE"):
        g, o = run(ARM_DEF[arm], U); s = L.summ(o); out[arm] = s
        if v1:
            ref = v1["arms"][arm]
            out[arm]["v1_median"] = ref["median"]; out[arm]["v1_p_success"] = ref["p_success"]
            out[arm]["z_success"] = (s["p_success"] - ref["p_success"]) / np.sqrt(s["se"] ** 2 + ref["se_success"] ** 2)
    tick("equiv done"); return out


def arms_all():
    d = {k: dict(v) for k, v in ARM_DEF.items()}
    P1 = ARM_DEF["P1_dual_lock"]; R1 = ARM_DEF["R1_reference_PE"]
    d["M0_matrix_only"] = dict(P1, pe=0.0, inh=0.0)
    d["R1_half_dose"] = dict(R1, pe=0.5)
    d["P1_full_dose"] = dict(P1, pe=1.0)
    for k, v in (("delta", 1.0), ("k_dil", 1.5), ("mism", 0.5), ("k_rel", 2.0)):
        d["P1_minus_%s" % k] = dict(P1, **{k: v})
    d["R1_plus_confinement_bundle"] = dict(R1, delta=0.30, k_dil=0.8, mism=0.15, k_rel=1.2)
    for k, v in (("delta", 0.30), ("k_dil", 0.8), ("mism", 0.15), ("k_rel", 1.2), ("prot", 0.70), ("k_dep", 3.0)):
        d["R1_plus_%s" % k] = dict(R1, **{k: v})
    d["R1_plus_halfPE_efflux"] = dict(R1, pe=0.5, inh=0.033)
    return d


def stage_paired():
    U = L.uniforms(N_MAIN, S0 + 2)
    arms = arms_all(); res = {"arms": {}, "paired_vs_R1": {}, "paired_vs_P1": {}}; store = {}
    for k, spec in arms.items():
        g, o = run(spec, U); store[k] = o; res["arms"][k] = L.summ(o); res["arms"][k]["closure_err"] = float(np.max(np.abs(o["closure"] - 1)))
    ok = {k: L.success(o) for k, o in store.items()}
    for k in arms:
        res["paired_vs_R1"][k] = dict(d_success=paired_ci(ok[k], ok["R1_reference_PE"]))
        res["paired_vs_P1"][k] = dict(d_success=paired_ci(ok[k], ok["P1_dual_lock"]))
    a, b = store["P1_dual_lock"], store["R1_reference_PE"]
    res["dominance_P1_vs_R1"] = dict(
        p_F_higher=float(np.mean(a["F"] > b["F"])), p_DH_lower=float(np.mean(a["DH"] < b["DH"])),
        p_both=float(np.mean((a["F"] > b["F"]) & (a["DH"] < b["DH"]))),
        p_P1_succeeds_R1_fails=float(np.mean(ok["P1_dual_lock"] & ~ok["R1_reference_PE"])),
        p_R1_succeeds_P1_fails=float(np.mean(~ok["P1_dual_lock"] & ok["R1_reference_PE"])),
        median_ratio_F=float(np.median(a["F"] / np.maximum(b["F"], 1e-12))), median_DH_ratio=float(np.median(a["DH"] / np.maximum(b["DH"], 1e-12))))
    res["samples"] = {k: dict(F=store[k]["F"][:3000].tolist(), DH=store[k]["DH"][:3000].tolist(), trev=(store[k]["t_last"][:3000]).tolist())
                      for k in ("R1_reference_PE", "P1_dual_lock", "P5_no_self_limiting_PE")}
    tick("paired done"); return res


def stage_converge():
    out = {}
    U = L.uniforms(100000, S0 + 3)
    for arm in ("R1_reference_PE", "P1_dual_lock"):
        g, o = run(ARM_DEF[arm], U); F = o["F"]; ok = L.success(o)
        cps = [250, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
        out["running_" + arm] = [dict(n=c, median_F=float(np.median(F[:c])), p_success=float(ok[:c].mean()), p_F_ge_4=float(np.mean(F[:c] >= M.F_TARGET))) for c in cps]
        nb = 20; bm = ok[:100000].reshape(nb, -1).mean(1)
        out["batch_means_" + arm] = dict(mean=float(bm.mean()), sd=float(bm.std(ddof=1)), se=float(bm.std(ddof=1) / np.sqrt(nb)))
    reps = {}
    for arm in ("R1_reference_PE", "P1_dual_lock"):
        med, ps = [], []
        for r in range(20):
            g, o = run(ARM_DEF[arm], L.uniforms(N_MAIN, S0 + 100 + r)); med.append(np.median(o["F"])); ps.append(L.success(o).mean())
        reps[arm] = dict(median_F_mean=float(np.mean(med)), median_F_sd=float(np.std(med, ddof=1)), p_success_mean=float(np.mean(ps)), p_success_sd=float(np.std(ps, ddof=1)),
                         median_F_range=[float(min(med)), float(max(med))], p_success_range=[float(min(ps)), float(max(ps))], n_reps=20, n_per_rep=N_MAIN)
    out["replicates"] = reps
    U5 = L.uniforms(5000, S0 + 4); out["dt"] = {}
    for dt in (0.04, 0.02, 0.01, 0.005, 0.0025):
        g, o = run(ARM_DEF["P1_dual_lock"], U5, dt=dt); out["dt"][str(dt)] = dict(median_F=float(np.median(o["F"])), p_success=float(L.success(o).mean()))
    out["horizon"] = {}
    for te in (4.0, 6.0, 8.0, 12.0, 24.0):
        for arm in ("R1_reference_PE", "P1_dual_lock"):
            g, o = run(ARM_DEF[arm], U5, t_end=te); out["horizon"].setdefault(arm, {})[str(te)] = dict(median_F=float(np.median(o["F"])), p_success=float(L.success(o).mean()))
    tick("converge done"); return out


def stage_corr():
    out = {"scenarios": {}, "sweep_rho_EC50_Ctox": {}}
    def both(corr, tag):
        U = L.uniforms(N_MAIN, S0 + 5, corr)
        r = {}
        for arm in ("R1_reference_PE", "P1_dual_lock", "P5_no_self_limiting_PE"):
            g, o = run(ARM_DEF[arm], U); r[arm] = L.summ(o)
        return r
    out["scenarios"]["independent"] = both(None, "ind")
    out["scenarios"]["base_correlated"] = both(CORR_BASE, "base")
    out["scenarios"]["base_correlated"]["definition"] = {"%s|%s" % k: v for k, v in CORR_BASE.items()}
    for rho in (0.0, 0.3, 0.6, 0.8, 0.9):
        out["sweep_rho_EC50_Ctox"][str(rho)] = both({("EC50", "C_tox"): rho} if rho else None, "sweep")
    out["sweep_other"] = {}
    for tag, cr in (("rho(C_tox,t_rev)=-0.8", {("C_tox", "t_rev"): -0.8}), ("rho(EC50,Emax)=-0.8", {("EC50", "Emax"): -0.8}),
                    ("adverse bundle: rho(EC50,C_tox)=0.8, rho(EC50,Emax)=-0.5, rho(C_tox,t_rev)=-0.6", {("EC50", "C_tox"): 0.8, ("EC50", "Emax"): -0.5, ("C_tox", "t_rev"): -0.6}),
                    ("favourable bundle: rho(EC50,C_tox)=-0.5, rho(EC50,Emax)=0.3, rho(C_tox,t_rev)=0.4", {("EC50", "C_tox"): -0.5, ("EC50", "Emax"): 0.3, ("C_tox", "t_rev"): 0.4})):
        out["sweep_other"][tag] = both(cr, tag)
    # realised rank correlations
    U = L.uniforms(50000, S0 + 6, CORR_BASE)
    out["realised_spearman"] = {"%s|%s" % k: float(spearmanr(U[:, L.IDX[k[0]]], U[:, L.IDX[k[1]]]).statistic) for k in CORR_BASE}
    tick("corr done"); return out


def _sobol_fun(spec, p0v):
    def f(U):
        g, o = L.run_u(U, spec, p0v)
        return {"log10_F": np.log10(np.maximum(o["F"], 1e-6)), "success": L.success(o).astype(float), "log10_DH": np.log10(np.maximum(o["DH"], 1e-6))}
    return f


def stage_sobol():
    out = {}
    for arm in ("P1_dual_lock", "R1_reference_PE"):
        out[arm] = L.sobol(_sobol_fun(ARM_DEF[arm], p0()), 13, S0 + 7, groups=L.GROUPS, nboot=200)
        tick("sobol %s done" % arm)
    out["n_base"] = 2 ** 13; out["n_evals_per_arm"] = 2 ** 13 * (len(L.DIMS) + 2 + len(L.GROUPS))
    return out


def stage_pareto():
    doses = [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    U = L.uniforms(10000, S0 + 8)
    fam = {"R1 reference (SNAC-like)": ARM_DEF["R1_reference_PE"], "R2 reference + efflux module": ARM_DEF["R2_reference_PE_efflux"], "P1 Dual-Lock": ARM_DEF["P1_dual_lock"],
           "P3 Dual-Lock without co-localisation": ARM_DEF["P3_no_colocalisation"], "P5 Dual-Lock without self-limiting PE": ARM_DEF["P5_no_self_limiting_PE"]}
    out = {"doses": doses, "families": {}}
    for name, spec in fam.items():
        out["families"][name] = []
        for d in doses:
            g, o = run(spec, U, pe=d); s = L.summ(o); s["dose"] = d; s["p_unsafe"] = 1 - s["p_safe"]; out["families"][name].append(s)
    tick("pareto done"); return out


def stage_fals():
    deltas = [0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.20]
    doses = [0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00, 1.50, 2.00]
    kdils = [0.2, 0.4, 0.8, 1.5, 2.5, 4.0]
    U = L.uniforms(5000, S0 + 9)
    P1 = ARM_DEF["P1_dual_lock"]; res = {"deltas": deltas, "doses": doses, "k_dils": kdils}
    keys = ("p_success", "median", "p_safe", "median_DH", "p_F_ge_4")
    res["grid_delta_dose"] = {k: np.zeros((len(doses), len(deltas))).tolist() for k in keys}
    G = {k: np.zeros((len(doses), len(deltas))) for k in keys}
    for j, dl in enumerate(deltas):
        for i, d in enumerate(doses):
            g, o = run(dict(P1, delta=dl), U, pe=d); s = L.summ(o)
            for k in keys:
                G[k][i, j] = s[k]
    res["grid_delta_dose"] = {k: v.tolist() for k, v in G.items()}
    H = {k: np.zeros((len(kdils), len(deltas))) for k in keys}
    for j, dl in enumerate(deltas):
        for i, kd in enumerate(kdils):
            g, o = run(dict(P1, delta=dl, k_dil=kd), U); s = L.summ(o)
            for k in keys:
                H[k][i, j] = s[k]
    res["grid_delta_kdil"] = {k: v.tolist() for k, v in H.items()}
    # predefined benefit region B: median F >= 4 % and P(safe) >= 0.80 ; failure of the architecture at a given delta = no dose reaches B
    inB = (G["median"] >= M.F_TARGET) & (G["p_safe"] >= 0.80)
    res["region_B"] = inB.tolist()
    res["max_p_success_by_delta"] = G["p_success"].max(0).tolist()
    res["best_dose_by_delta"] = [doses[int(i)] for i in G["p_success"].argmax(0)]
    res["delta_max_for_B"] = float(max([dl for j, dl in enumerate(deltas) if inB[:, j].any()] or [0.0]))
    res["delta_max_for_p30"] = float(max([dl for j, dl in enumerate(deltas) if G["p_success"][:, j].max() >= 0.30] or [0.0]))
    res["definitions"] = dict(B="median F >= 4 % and P(safe) >= 0.80 for some dose", p30="max over dose of P(success) >= 0.30")
    tick("fals done"); return res


def stage_struct():
    out = {"variants": {}}
    n = 2000; U = L.uniforms(n, S0 + 10); UD = L.uniforms(n, S0 + 11)[:, :3]; UL = L.uniforms(n, S0 + 12)[:, :2]
    specs = {k: ARM_DEF[k] for k in ("R1_reference_PE", "P1_dual_lock", "P2_no_efflux_module", "P3_no_colocalisation", "P5_no_self_limiting_PE")}
    p0_l = p0()

    def summarize(runner, p0v):
        r = {}
        for k, sp in specs.items():
            o = runner(sp, p0v); r[k] = L.summ(o)
        r["R1_over_P1_ratio_median_F"] = r["P1_dual_lock"]["median"] / r["R1_reference_PE"]["median"]
        r["delta_success_P1_minus_R1"] = r["P1_dual_lock"]["p_success"] - r["R1_reference_PE"]["p_success"]
        return r

    # (a) lumped, recalibrated sampler on the same U
    lum = lambda sp, p: L.run_u(U, sp, p)[1]
    pc = L.calibrate_generic(lambda p: lum(specs["R1_reference_PE"], p), iters=24, lo=1e-4, hi=5.0)
    out["variants"]["lumped (base)"] = dict(p0=pc, p0_over_base=pc / p0_l, **summarize(lum, pc))
    # (b) lumped + stagnant layer, recalibrated
    lay = lambda sp, p: L.run_u(U, sp, p, layer=UL)[1]
    pc = L.calibrate_generic(lambda p: lay(specs["R1_reference_PE"], p), iters=24, lo=1e-4, hi=5.0)
    out["variants"]["lumped + stagnant layer"] = dict(p0=pc, p0_over_base=pc / p0_l, **summarize(lay, pc))
    # layer without recalibration (P0 kept at the base value): shows the cap on the enhanced permeability
    out["variants"]["lumped + stagnant layer, P0 not recalibrated"] = dict(p0=p0_l, p0_over_base=1.0, **summarize(lay, p0_l))
    tick("struct lumped variants done")
    # (c) spatial model, recalibrated to the same 1.2 % anchor, for several effective dispersion multipliers
    for mix in (3, 10, 30, 100, 1000, 1e6):
        pdeR = lambda sp, p, mix=mix: L.run_pde(U, sp, p, UD, mix=mix)[1]
        fR1 = lambda p, mix=mix: L.run_pde(U[:1000], specs["R1_reference_PE"], p, UD[:1000], mix=mix)[1]
        hi = 40.0
        top = np.median(fR1(hi)["F"])
        if top < M.F_BASELINE_MID:
            out["variants"]["spatial, dispersion x%g" % mix] = dict(calibratable=False, max_median_F_R1_at_P0_40=float(top), note="reference arm cannot reach 1.2 % even at P0 = 40 cm/h")
            tick("spatial mix %g: not calibratable" % mix); continue
        pc = L.calibrate_generic(fR1, iters=20, lo=1e-4, hi=hi)
        r = summarize(pdeR, pc)
        out["variants"]["spatial, dispersion x%g" % mix] = dict(calibratable=True, p0=pc, p0_over_base=pc / p0_l, **r)
        # same variant without recalibration: P0 fixed at the base value
        r2 = summarize(pdeR, p0_l)
        out["variants"]["spatial, dispersion x%g, P0 not recalibrated" % mix] = dict(calibratable=True, p0=p0_l, p0_over_base=1.0, **r2)
        tick("spatial mix %g done (P0/base = %.2f)" % (mix, pc / p0_l))
    # diffusion time scales for the reader
    out["diffusion_times_h"] = dict(P1_peptide=(0.3 ** 2) / (L.DIFF["P_cm2_s"] * 3600), R1_peptide=(1.0 ** 2) / (L.DIFF["P_cm2_s"] * 3600),
                                    layer_peptide=(L.LAYER["h_s_cm"] ** 2) / (L.LAYER["D_cm2_s"] * 3600))
    out["n_draws"] = n
    return out


def stage_abc():
    n = 100000; U = L.uniforms(n, S0 + 13)
    g1, oR = run(ARM_DEF["R1_reference_PE"], U); g2, oP = run(ARM_DEF["P1_dual_lock"], U)
    out = {"n": n}
    names = ["P0", "k_deg", "S_post", "ER0", "A_m", "EC50", "Emax", "h", "C_tox", "t_rev", "IC50"]
    for tag, (lo, hi) in {"prior (no conditioning)": (0, 1), "R1 F in [0.4, 2.0] %": (0.004, 0.02), "R1 F in [0.4, 1.0] % (oral-semaglutide-like)": (0.004, 0.01), "R1 F in [0.8, 1.6] %": (0.008, 0.016)}.items():
        acc = (oR["F"] >= lo) & (oR["F"] <= hi)
        r = dict(n_accepted=int(acc.sum()), frac=float(acc.mean()),
                 P1=dict(median=float(np.median(oP["F"][acc])), p_success=float(L.success(oP)[acc].mean()), p_F_ge_4=float(np.mean(oP["F"][acc] >= M.F_TARGET)),
                         p05=float(np.percentile(oP["F"][acc], 5)), p95=float(np.percentile(oP["F"][acc], 95))),
                 R1=dict(median=float(np.median(oR["F"][acc])), p_success=float(L.success(oR)[acc].mean())),
                 shift={k: float(np.median(g1[k][acc]) / np.median(g1[k])) for k in names})
        out[tag] = r
    tick("abc done"); return out


def stage_valid():
    out = {}
    U = L.uniforms(N_MAIN, S0 + 14)
    g, oR = run(ARM_DEF["R1_reference_PE"], U); g2, oP = run(ARM_DEF["P1_dual_lock"], U)
    lnF = np.log(np.maximum(oR["F"], 1e-12))
    out["R1"] = dict(median=float(np.median(oR["F"])), sd_lnF=float(lnF.std()), geometric_cv=float(np.sqrt(np.exp(lnF.var()) - 1)),
                     frac_in_0p4_2p0=float(np.mean((oR["F"] >= 0.004) & (oR["F"] <= 0.02))), frac_in_0p4_1p0=float(np.mean((oR["F"] >= 0.004) & (oR["F"] <= 0.01))),
                     frac_gt_6=float(np.mean(oR["F"] > 0.06)), frac_gt_10=float(np.mean(oR["F"] > 0.10)), frac_gt_20=float(np.mean(oR["F"] > 0.20)))
    out["P1"] = dict(median=float(np.median(oP["F"])), frac_gt_6=float(np.mean(oP["F"] > 0.06)), frac_gt_10=float(np.mean(oP["F"] > 0.10)),
                     frac_gt_20=float(np.mean(oP["F"] > 0.20)), frac_gt_30=float(np.mean(oP["F"] > 0.30)))
    out["observed"] = dict(oral_semaglutide_BAV_pct=0.8, oral_semaglutide_within_subject_CV_pct=137.0, typical_enhancer_BAV_pct=1.0, MEDI7219_dog_BAV_pct=6.0)
    out["observed_equiv_sd_ln"] = float(np.sqrt(np.log(1 + 1.37 ** 2)))
    # monotonic / qualitative predictions
    q = {}
    for tag, mult in (("k_dil x0.5", 0.5), ("k_dil x1", 1.0), ("k_dil x2", 2.0)):
        gg = L.inputs_from_u(U[:8000], ARM_DEF["R1_reference_PE"], p0()); gg["k_dil"] = gg["k_dil"] * mult
        q["R1 " + tag] = float(np.median(M.integrate(gg, "R1_reference_PE")["F"]))
    for dl in (0.5, 1.0, 2.0):
        gg = L.inputs_from_u(U[:8000], dict(ARM_DEF["R1_reference_PE"], delta=dl), p0()); q["R1 delta=%g cm" % dl] = float(np.median(M.integrate(gg, "R1_reference_PE")["F"]))
    for d in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        gg = L.inputs_from_u(U[:8000], ARM_DEF["R1_reference_PE"], p0()); q["R1 PE dose %g mmol" % d] = float(np.median(M.integrate(gg, "R1_reference_PE", pe=d)["F"]))
    out["qualitative"] = q
    # transience: fraction of final uptake already absorbed by time t (median over draws)
    tr = {}
    for arm in ("R1_reference_PE", "P1_dual_lock"):
        fin = L.run_u(U[:5000], ARM_DEF[arm], p0())[1]["Abs"]
        tr[arm] = {}
        for t in (0.5, 1.0, 2.0, 4.0):
            a_t = L.run_u(U[:5000], ARM_DEF[arm], p0(), t_end=t)[1]["Abs"]
            tr[arm][str(t)] = float(np.median(a_t / np.maximum(fin, 1e-15)))
    out["fraction_absorbed_by_t"] = tr
    tick("valid done"); return out



def _models(U, UL, US):
    """The structural / constraint variants compared in the revision. Each returns (g, outputs) for a spec and baseline permeability p0v."""
    return {
        "lumped (original)": lambda sp, p, **kw: L.run_u(U, sp, p, **kw),
        "uptake-constrained": lambda sp, p, **kw: L.run_u(U, sp, p, layer=UL, sol=US, **kw),
    }


def stage_constrained():
    """Matrix-only artefact and biologically constrained recalibration: stagnant layer (serial resistance) + peptide-solubility limit on the dissolved amount,
    each variant recalibrated to the same 1.2 % anchor; matrix-only gain, volume-only effects, upper tail, damage-at-efficacy and tail-rejection analyses."""
    n = 20000; U = L.uniforms(n, S0 + 20); UL = L.uniforms(n, S0 + 21)[:, :2]; US = L.uniforms(n, S0 + 22)[:, :2]
    P1, R1 = ARM_DEF["P1_dual_lock"], ARM_DEF["R1_reference_PE"]
    specs = {"R0_peptide_alone": ARM_DEF["R0_peptide_alone"], "M0_matrix_only": dict(P1, pe=0.0, inh=0.0), "R1_reference_PE": R1, "P1_dual_lock": P1,
             "R1_plus_delta": dict(R1, delta=0.30), "P1_minus_delta": dict(P1, delta=1.0), "P3_no_colocalisation": ARM_DEF["P3_no_colocalisation"], "P5_no_self_limiting_PE": ARM_DEF["P5_no_self_limiting_PE"],
             "R1_plus_confinement_bundle": dict(R1, delta=0.30, k_dil=0.8, mism=0.15, k_rel=1.2)}
    variants = {"lumped (original)": dict(), "stagnant layer": dict(layer=UL), "solubility limit": dict(sol=US), "uptake-constrained (layer + solubility)": dict(layer=UL, sol=US),
                "uptake-constrained, damage at efficacious concentration": dict(layer=UL, sol=US, mods={"C_tox": 25.0 / 80.0})}
    out = {"variants": {}, "n": n}
    for name, kw in variants.items():
        calkw = {k: v for k, v in kw.items() if k != "mods"}
        fcal = lambda p, kw=kw: L.run_u(U[:8000], R1, p, **{k: (v[:8000] if hasattr(v, "shape") else v) for k, v in kw.items()})[1]
        pc_ = L.calibrate_generic(fcal, iters=24, lo=1e-4, hi=5.0)
        res = {}; store = {}
        for k, sp in specs.items():
            g, o = L.run_u(U, sp, pc_, **kw); store[k] = (g, o); res[k] = L.summ(o); res[k]["frac_gt_10"] = float(np.mean(o["F"] > 0.10)); res[k]["frac_gt_6"] = float(np.mean(o["F"] > 0.06))
        res["p0"] = pc_; res["p0_over_base"] = pc_ / p0()
        res["matrix_only_gain"] = res["M0_matrix_only"]["median"] / res["R0_peptide_alone"]["median"]
        # observation-conditioned and tail-rejection analyses within this variant
        oR, oP = store["R1_reference_PE"][1], store["P1_dual_lock"][1]
        okP = L.success(oP)
        acc = (oR["F"] >= 0.004) & (oR["F"] <= 0.02)
        accs = (oR["F"] >= 0.004) & (oR["F"] <= 0.01)
        cred = (oR["F"] <= 0.10) & (oP["F"] <= 0.10)
        res["abc_0.4_2.0"] = dict(frac=float(acc.mean()), P1_success=float(okP[acc].mean()), P1_median=float(np.median(oP["F"][acc])), P1_frac_gt_10=float(np.mean(oP["F"][acc] > 0.10)))
        res["abc_0.4_1.0"] = dict(frac=float(accs.mean()), P1_success=float(okP[accs].mean()), P1_median=float(np.median(oP["F"][accs])))
        res["tail_rejected_F_le_10"] = dict(frac=float(cred.mean()), P1_success=float(okP[cred].mean()), R1_success=float(L.success(oR)[cred].mean()), P1_median=float(np.median(oP["F"][cred])))
        both = acc & cred
        res["abc_and_tail"] = dict(frac=float(both.mean()), P1_success=float(okP[both].mean()), P1_median=float(np.median(oP["F"][both])))
        out["variants"][name] = res
        tick("constrained: %s done (P0 x%.2f)" % (name, pc_ / p0()))
    return out


def stage_revers():
    """Reversibility endpoints: sensitivity of success to the damage-hours limit, the recovery limit and the potency/toxicity ratio C_tox/EC50."""
    n = 20000; U = L.uniforms(n, S0 + 30); res = {"scenarios": {}, "ctox_ratio": {}}
    arms = {"R1_reference_PE": ARM_DEF["R1_reference_PE"], "P1_dual_lock": ARM_DEF["P1_dual_lock"], "P5_no_self_limiting_PE": ARM_DEF["P5_no_self_limiting_PE"]}
    for mult in (1.0, 0.3, 0.1, 0.03):
        res["ctox_ratio"][str(mult)] = {"C_tox_over_EC50_median": 80.0 * mult / 25.0}
        for k, sp in arms.items():
            g, o = L.run_u(U, sp, p0(), mods={"C_tox": mult})
            F = o["F"]; DH = o["DH"]; trev = o["t_last"] + g["t_rev"]
            d = dict(median_DH=float(np.median(DH)), p90_DH=float(np.percentile(DH, 90)), p_DH_gt_0p5=float(np.mean(DH > 0.5)), p_success=float(L.success(o).mean()), p_safe=float(o["safe"].mean()))
            d["grid"] = {"%s|%s" % (dh, tr): float(np.mean((F >= M.F_TARGET) & (DH <= dh) & (trev <= tr))) for dh in (0.5, 0.1, 0.02, 0.005, 0.001) for tr in (2.0, 3.0, 4.0, 6.0)}
            d["p_recov_le"] = {str(tr): float(np.mean(trev <= tr)) for tr in (1.0, 2.0, 3.0, 4.0, 6.0)}
            res["ctox_ratio"][str(mult)][k] = d
    tick("revers done"); return res


def stage_optim():
    """Out-of-sample evaluation of a formulation optimized on training draws: design variables k_rel, enhancer dose and wash-out are chosen on a training set
    and evaluated, together with the nominal design, on an independent test set (new seed), for the original lumped model and the uptake-constrained model."""
    P1 = ARM_DEF["P1_dual_lock"]; out = {}
    krels = [1.2, 1.6, 2.0, 2.5, 3.0, 4.0]; doses = [0.35, 0.5, 0.75, 1.0]; kdils = [0.4, 0.8, 1.5]
    for model, kw in (("lumped (original)", "lumped"), ("uptake-constrained", "constrained")):
        Utr = L.uniforms(8000, S0 + 40); Ute = L.uniforms(20000, S0 + 41)
        kws = lambda U_, off: ({} if kw == "lumped" else dict(layer=L.uniforms(len(U_), S0 + off)[:, :2], sol=L.uniforms(len(U_), S0 + off + 1)[:, :2]))
        ktr, kte = kws(Utr, 50), kws(Ute, 52)
        pc_ = p0()
        if kw == "constrained":
            Uc = L.uniforms(8000, S0 + 20); pc_ = L.calibrate_generic(lambda p: L.run_u(Uc, ARM_DEF["R1_reference_PE"], p, layer=L.uniforms(8000, S0 + 21)[:, :2], sol=L.uniforms(8000, S0 + 22)[:, :2])[1], iters=24, lo=1e-4, hi=5.0)
        grid = []
        for kr in krels:
            for d_ in doses:
                for kd in kdils:
                    sp = dict(P1, k_rel=kr, k_dil=kd); g, o = L.run_u(Utr, sp, pc_, pe=d_, **ktr); grid.append((float(L.success(o).mean()), kr, d_, kd, float(np.mean(o["F"] > 0.10)), float(o["safe"].mean())))
        grid.sort(reverse=True); best = grid[0]
        def ev(sp, pe_, U_, kk):
            g, o = L.run_u(U_, sp, pc_, pe=pe_, **kk); s = L.summ(o); s["frac_gt_10"] = float(np.mean(o["F"] > 0.10)); return s
        nom_tr = ev(P1, 0.5, Utr, ktr); nom_te = ev(P1, 0.5, Ute, kte)
        opt_sp = dict(P1, k_rel=best[1], k_dil=best[3])
        opt_tr = ev(opt_sp, best[2], Utr, ktr); opt_te = ev(opt_sp, best[2], Ute, kte)
        ref_te = ev(ARM_DEF["R1_reference_PE"], ARM_DEF["R1_reference_PE"]["pe"], Ute, kte)
        # paired difference on the test set
        g1, o1 = L.run_u(Ute, P1, pc_, pe=0.5, **kte); g2, o2 = L.run_u(Ute, opt_sp, pc_, pe=best[2], **kte)
        cred = [t for t in grid if t[4] <= 0.30 and t[1] <= 2.5]
        bc = cred[0]; cred_sp = dict(P1, k_rel=bc[1], k_dil=bc[3])
        cred_tr = ev(cred_sp, bc[2], Utr, ktr); cred_te = ev(cred_sp, bc[2], Ute, kte)
        g3, o3 = L.run_u(Ute, cred_sp, pc_, pe=bc[2], **kte)
        out[model + "|credible"] = dict(best=dict(k_rel=bc[1], dose=bc[2], k_dil=bc[3], train_success=bc[0]), optimized_train=cred_tr, optimized_test=cred_te, nominal_test=nom_te, nominal_train=nom_tr, reference_test=ref_te,
                                        paired_opt_minus_nominal_test=paired_ci(L.success(o3), L.success(o1)), rule="maximize training P(success) subject to P(F>10%) <= 0.30 and k_rel <= 2.5 per h")
        out[model] = dict(best=dict(k_rel=best[1], dose=best[2], k_dil=best[3], train_success=best[0]), top5=[dict(train_success=a, k_rel=b, dose=c, k_dil=d) for a, b, c, d, e_, f_ in grid[:5]],
                          nominal_train=nom_tr, nominal_test=nom_te, optimized_train=opt_tr, optimized_test=opt_te, reference_test=ref_te,
                          paired_opt_minus_nominal_test=paired_ci(L.success(o2), L.success(o1)), n_train=8000, n_test=20000, n_grid=len(grid))
        tick("optim %s done: best %s" % (model, best))
    return out



def stage_decomp():
    """Decomposition of the matrix-only gain: which formulation parameter of the matrix (niche thickness, wash-out, protection, release) produces the gain over peptide alone,
    in the original lumped model and in the uptake-constrained model (same recalibrated baselines as stage 'constrained')."""
    cons = json.load(open(os.path.join(OUT, "ext_constrained.json")))["variants"]
    n = 20000; U = L.uniforms(n, S0 + 20); UL = L.uniforms(n, S0 + 21)[:, :2]; US = L.uniforms(n, S0 + 22)[:, :2]
    R0 = ARM_DEF["R0_peptide_alone"]; P1 = ARM_DEF["P1_dual_lock"]
    # matrix parameters applied to peptide alone (no enhancer, no modulator); the protection of R0 is 0, that of the matrix 0.70
    def arm(**kw):
        return dict(R0, pe=0.0, inh=0.0, **kw)
    arms = {"R0 peptide alone": arm(), "R0 + thin niche (delta 0.30)": arm(delta=0.30), "R0 + slow wash-out (0.8)": arm(k_dil=0.8), "R0 + protection 0.70": arm(prot=0.70),
            "R0 + niche + wash-out": arm(delta=0.30, k_dil=0.8), "R0 + niche + wash-out + synchronous release": arm(delta=0.30, k_dil=0.8, mism=0.15, k_rel=1.2),
            "M0 matrix only (all matrix parameters)": dict(P1, pe=0.0, inh=0.0), "M0 without protection": dict(P1, pe=0.0, inh=0.0, prot=0.0), "M0 without niche thinning": dict(P1, pe=0.0, inh=0.0, delta=1.0)}
    out = {}
    for name, kw in (("lumped (original)", {}), ("uptake-constrained (layer + solubility)", dict(layer=UL, sol=US))):
        p0v = cons[name]["p0"]; res = {}
        for k, sp in arms.items():
            g, o = L.run_u(U, sp, p0v, **kw); res[k] = dict(median=float(np.median(o["F"])), p95=float(np.percentile(o["F"], 95)))
        base = res["R0 peptide alone"]["median"]
        for k in res:
            res[k]["gain"] = res[k]["median"] / base
        out[name] = res
    tick("decomp done"); return out


STAGES = ["equiv", "paired", "converge", "corr", "sobol", "pareto", "fals", "struct", "abc", "valid", "constrained", "revers", "optim", "decomp"]


# ------------------------------------------------------------------------------------------------ verification tests
def self_tests():
    ok = True

    def chk(name, cond, info=""):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + name + (" [%s]" % info if info else ""), flush=True)
        ok = ok and bool(cond)

    # 1 equivalence of the U-sampler with the v1 sampler (same marginals)
    pz = p0(); U = L.uniforms(20000, 5)
    g, o = run(ARM_DEF["P1_dual_lock"], U); s = L.summ(o)
    _, o1 = M.run_arm(20000, "P1_dual_lock", pz, seed_offset=100 + 3); s1 = M.summarise(o1)
    z = abs(s["p_success"] - s1["p_success"]) / np.sqrt(s["se"] ** 2 + s1["se_success"] ** 2)
    chk("U-sampler reproduces the v1 P1 success probability (|z| < 3)", z < 3, "%.3f vs %.3f, z=%.2f" % (s["p_success"], s1["p_success"], z))
    chk("U-sampler median F within 6 % of v1", abs(s["median"] / s1["median"] - 1) < 0.06, "%.4f vs %.4f" % (s["median"], s1["median"]))
    # 2 mass balance
    chk("lumped mass balance closes (1e-9)", np.max(np.abs(o["closure"] - 1)) < 1e-9)
    # 3 Sobol estimator on the Ishigami function
    def ish(U_):
        x = -np.pi + 2 * np.pi * U_[:, :3]
        y = np.sin(x[:, 0]) + 7 * np.sin(x[:, 1]) ** 2 + 0.1 * x[:, 2] ** 4 * np.sin(x[:, 0])
        return {"y": y}
    r = L.sobol(ish, 14, 3, groups=None, nboot=10)["y"]["params"]
    exact = {"P0": (0.3139, 0.5576), "k_deg": (0.4424, 0.4424), "S_post": (0.0, 0.2437)}
    err = max(max(abs(r[k]["S1"] - v[0]), abs(r[k]["ST"] - v[1])) for k, v in exact.items())
    chk("Sobol estimators recover the Ishigami indices (max error < 0.04)", err < 0.04, "max error %.3f" % err)
    chk("Sobol: inert inputs have indices ~ 0", max(abs(r[k]["ST"]) for k in L.DIMS[3:]) < 0.03)
    # 4 Gaussian copula realises the requested rank correlation
    Uc = L.uniforms(60000, 8, {("EC50", "C_tox"): 0.6})
    rho_s = spearmanr(Uc[:, L.IDX["EC50"]], Uc[:, L.IDX["C_tox"]]).statistic
    chk("copula Spearman matches (6/pi) asin(rho/2)", abs(rho_s - 6 / np.pi * np.arcsin(0.3)) < 0.01, "%.3f" % rho_s)
    # 5 spatial model: mixing limit and mass balance
    n = 1200; Us = L.uniforms(n, 9); UD = L.uniforms(n, 10)[:, :3]
    gl, ol = L.run_u(Us, ARM_DEF["P1_dual_lock"], pz)
    gp, op = L.run_pde(Us, ARM_DEF["P1_dual_lock"], pz, UD, mix=1e6)
    chk("spatial model -> lumped model in the well-mixed limit (median F within 3 %)", abs(np.median(op["F"]) / np.median(ol["F"]) - 1) < 0.03, "%.4f vs %.4f" % (np.median(op["F"]), np.median(ol["F"])))
    gp2, op2 = L.run_pde(Us, ARM_DEF["P1_dual_lock"], pz, UD, mix=50)
    chk("spatial model peptide mass balance closes (1e-7)", np.max(np.abs(op2["closure"] - 1)) < 1e-7, "%.1e" % np.max(np.abs(op2["closure"] - 1)))
    gp3, op3 = L.run_pde(Us, ARM_DEF["P1_dual_lock"], pz, UD, mix=50, nx=48, dt=0.01)
    chk("spatial model mesh/time convergence (median F within 3 %)", abs(np.median(op3["F"]) / np.median(op2["F"]) - 1) < 0.03, "%.4f vs %.4f" % (np.median(op3["F"]), np.median(op2["F"])))
    chk("spatial model: slower mixing lowers F (unrecalibrated, same draws)", np.median(L.run_pde(Us, ARM_DEF["P1_dual_lock"], pz, UD, mix=5)[1]["F"]) < np.median(op2["F"]))
    # 6 stagnant layer
    UL = L.uniforms(n, 11)[:, :2]
    gL, oL = L.run_u(Us, ARM_DEF["P1_dual_lock"], pz, layer=UL)
    chk("stagnant layer never raises F (same draws)", np.all(oL["F"] <= ol["F"] + 1e-12))
    chk("lumped + layer mass balance closes (1e-9)", np.max(np.abs(oL["closure"] - 1)) < 1e-9)
    # 7 module logic with U inputs
    gx, ox = run(ARM_DEF["P1_dual_lock"], L.uniforms(3000, 12), p_sub=0.0)
    _, oy = run(dict(ARM_DEF["P1_dual_lock"], inh=0.0), L.uniforms(3000, 12), p_sub=0.0)
    chk("efflux module has no effect if the peptide is never a substrate", np.allclose(ox["F"], oy["F"]))
    _, oz = run(dict(ARM_DEF["P1_dual_lock"], pe=0.0, inh=0.0), L.uniforms(3000, 12))
    _, op1 = run(ARM_DEF["P1_dual_lock"], L.uniforms(3000, 12))
    chk("PE and modulator add a factor >= 4 to the matrix-only median F", np.median(op1["F"]) > 4 * np.median(oz["F"]), "%.4f vs %.4f" % (np.median(op1["F"]), np.median(oz["F"])))
    # 8 determinism
    a1 = run(ARM_DEF["P1_dual_lock"], L.uniforms(2000, 13))[1]["F"]; a2 = run(ARM_DEF["P1_dual_lock"], L.uniforms(2000, 13))[1]["F"]
    chk("deterministic with a fixed seed", np.array_equal(a1, a2))
    # 9 paired-difference helper
    d = paired_ci(np.ones(100) > 0, np.zeros(100) > 0)
    chk("paired_ci returns +1 for an all-true minus all-false contrast", abs(d[0] - 1) < 1e-12)
    # 10 solubility limit
    US2 = L.uniforms(1200, 15)[:, :2]
    g0, o0 = L.run_u(Us, ARM_DEF["P1_dual_lock"], pz)
    gS, oS = L.run_u(Us, ARM_DEF["P1_dual_lock"], pz, sol=US2)
    chk("solubility limit never raises F (same draws)", np.all(oS["F"] <= o0["F"] + 1e-12))
    chk("solubility limit: mass balance closes (1e-9)", np.max(np.abs(oS["closure"] - 1)) < 1e-9, "%.1e" % np.max(np.abs(oS["closure"] - 1)))
    gM, oM = L.run_u(Us, ARM_DEF["P1_dual_lock"], pz, mods={"C_tox": 0.1})
    chk("lowering C_tox raises damage-hours (same draws)", np.all(oM["DH"] >= o0["DH"] - 1e-12))
    print("ALL PASS" if ok else "SOME FAILED")
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--test" in args:
        sys.exit(0 if self_tests() else 1)
    todo = [a for a in args if not a.startswith("--")] or STAGES
    for st in todo:
        t = time.time(); res = globals()["stage_" + st](); save(st, res); tick("stage %s saved (%.0f s)" % (st, time.time() - t))
