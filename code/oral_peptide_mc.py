#!/usr/bin/env python3
"""
oral_peptide_mc.py  --  (v2.0: model unchanged from v1; tie-correct Spearman) Evidence-graded Monte Carlo screening model for oral
absolute bioavailability (F) of a 2-4 kDa peptide co-formulated with a
transcellular permeation enhancer (PE) and an efflux/transporter modulator.

Author : Leon Sandler, Independent Researcher
Style  : one file, numpy + matplotlib only, fixed seed, deterministic output,
         self-tests (--test), every number in the white paper read from results.json.

WHAT THIS IS
    A lumped, mechanistic prior-propagation model.  Parameter priors are
    author-assigned and evidence-graded (A-D, see PRIORS).  The model answers:
    "IF these priors hold, how likely is F >= 4 % with an acceptable and
    reversible epithelial exposure, which module contributes what, and what
    must an experiment show for the concept to survive?"

WHAT THIS IS NOT
    It is not experimental data.  It does not predict the bioavailability of
    any specific peptide.  Probabilities are conditional on the priors and on
    the model structure; they are NOT empirical probabilities of success.
    Calibration anchor: the reference arm is scaled so that its median F equals
    the Seeker's stated current baseline (1-2 %; mid-point 1.2 %).

MODEL (per Monte Carlo draw, explicit Euler, dt = 0.01 h, horizon 8 h)
    tablet pools : Pt (peptide), Et (PE), It (inhibitor)   -- first-order release
    micro-niche  : volume V = A_m * delta  (delta = niche thickness; matrix
                   confinement makes delta small)
        dP/dt = k_rel*Pt - (k_dil + k_deg*(1-Pr) + k_abs) * P
        dE/dt = k_rel*m*Et - (k_dil + k_dep) * E         (k_dep = enzymatic PE clearance)
        dI/dt = k_rel*m*It - k_dil * I
        Ce = E/V, Ci = I/V
        Hill  H = Ce^h / (EC50^h + Ce^h)      (Ce in mM = 1000*E[mmol]/V[mL])
        P_eff = P0 * (1 + Emax*H)                         (transcellular permeability)
        ER    = 1 + (ER0-1) / (1 + Ci/IC50)               (efflux ratio, peptide-substrate draws only)
        k_abs = P_eff / (delta * ER)
        F     = [integral k_abs*P dt] * S_post             (S_post = post-absorptive survival)
    exposure safety: damage-hours DH = integral Hill(Ce; C_tox, 3) dt
        SAFE  = (DH <= 0.5 h) and (t_clear + t_rev <= 4 h)
        t_clear = last time Ce > 0.1*EC50 ; t_rev = epithelial recovery time (prior)
    SUCCESS = (F >= 4 %) and SAFE      (challenge: >= 4-5 %, reversible within hours)

USAGE
    python oral_peptide_mc.py --test          # self-tests
    python oral_peptide_mc.py --run           # full analysis -> results/results.json + figures
    python oral_peptide_mc.py --run --n 5000  # quick run
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np

SEED = 20261031
N_DEFAULT = 20000
T_END = 8.0
DT = 0.01
F_TARGET = 0.04
F_STRETCH = 0.05
F_BASELINE_MID = 0.012      # Seeker baseline 1-2 %  -> calibration anchor (mid-point)
DH_MAX = 0.5                # h, tolerated damage-hours
T_REVERSE_MAX = 4.0         # h, challenge: "fully reversible within hours"

# --------------------------------------------------------------------------
# Evidence-graded priors.  grade: A measured, directly applicable literature value
# B measured for a close analogue   C order-of-magnitude literature range
# D author-assigned assumption (no direct source).  None is grade A here.
# kind: 'ln' lognormal(median, sigma_log)   'beta' (mean, sd)
# --------------------------------------------------------------------------
PRIORS = {
    "P0":       dict(kind="ln",   med=None, sig=0.60, unit="cm/h",  grade="C",
                     note="baseline apparent permeability of a hydrophilic 2-4 kDa peptide; median CALIBRATED to the Seeker baseline"),
    "k_deg":    dict(kind="ln",   med=2.0,  sig=0.50, unit="1/h",   grade="D",
                     note="luminal proteolysis rate in the niche without protection"),
    "S_post":   dict(kind="beta", mean=0.55, sd=0.15, unit="-",     grade="D",
                     note="survival after apical entry (intra-enterocyte handling, first pass)"),
    "ER0":      dict(kind="ln",   med=2.5,  sig=0.50, unit="-",     grade="C",
                     note="Caco-2 efflux ratio if the peptide is an efflux substrate"),
    "p_sub":    dict(kind="const", val=0.5, unit="-",               grade="D",
                     note="probability that the peptide is an efflux substrate (scenario axis)"),
    "A_m":      dict(kind="ln",   med=30.0, sig=0.30, unit="cm2",   grade="D",
                     note="mucosal contact area of the niche"),
    "EC50":     dict(kind="ln",   med=25.0, sig=0.50, unit="mM",    grade="C",
                     note="PE concentration for half-maximal permeability enhancement (millimolar PE activity, Maher 2016; Twarog 2019)"),
    "Emax":     dict(kind="ln",   med=80.0, sig=0.50, unit="-",     grade="D",
                     note="maximal fold-increase of peptide permeability by the PE"),
    "h":        dict(kind="ln",   med=2.0,  sig=0.25, unit="-",     grade="D",
                     note="Hill coefficient of PE effect"),
    "C_tox":    dict(kind="ln",   med=80.0, sig=0.40, unit="mM",    grade="D",
                     note="PE concentration for half-maximal epithelial damage"),
    "t_rev":    dict(kind="ln",   med=1.0,  sig=0.50, unit="h",     grade="D",
                     note="epithelial recovery time after PE falls below threshold"),
    "IC50":     dict(kind="ln",   med=0.5,  sig=0.80, unit="mM",    grade="C",
                     note="efflux-modulator potency (piperine 15-74 uM in Caco-2, Bhardwaj 2002; TPGS sub-mM, Collnot 2006)"),
}

# Design arms: medians of formulation-controlled parameters.
#   pe_mmol / inh_mmol : doses ; delta : niche thickness (cm) ; k_dil : niche washout (1/h)
#   mism : sigma_log of PE/peptide release mismatch ; prot : mean protection Pr (0-1)
#   k_dep: enzymatic PE clearance (1/h) ; k_rel : tablet release rate (1/h)
ARMS = {
    "R0_peptide_alone":      dict(pe=0.00, inh=0.000, delta=1.00, k_dil=1.50, mism=0.60, prot=0.00, k_dep=0.5, k_rel=2.0),
    "R1_reference_PE":       dict(pe=1.00, inh=0.000, delta=1.00, k_dil=1.50, mism=0.50, prot=0.50, k_dep=0.5, k_rel=2.0),
    "R2_reference_PE_efflux":dict(pe=1.00, inh=0.033, delta=1.00, k_dil=1.50, mism=0.50, prot=0.50, k_dep=0.5, k_rel=2.0),
    "P1_dual_lock":          dict(pe=0.50, inh=0.033, delta=0.30, k_dil=0.80, mism=0.15, prot=0.70, k_dep=3.0, k_rel=1.2),
    "P2_no_efflux_module":   dict(pe=0.50, inh=0.000, delta=0.30, k_dil=0.80, mism=0.15, prot=0.70, k_dep=3.0, k_rel=1.2),
    "P3_no_colocalisation":  dict(pe=0.50, inh=0.033, delta=1.00, k_dil=1.50, mism=0.50, prot=0.70, k_dep=3.0, k_rel=2.0),
    "P4_no_protection":      dict(pe=0.50, inh=0.033, delta=0.30, k_dil=0.80, mism=0.15, prot=0.50, k_dep=3.0, k_rel=1.2),
    "P5_no_self_limiting_PE":dict(pe=0.50, inh=0.033, delta=0.30, k_dil=0.80, mism=0.15, prot=0.70, k_dep=0.5, k_rel=1.2),
}


def _ln(rng, med, sig, n):
    return med * np.exp(sig * rng.standard_normal(n))


def _beta_ms(rng, mean, sd, n):
    var = sd ** 2
    k = mean * (1 - mean) / var - 1.0
    return rng.beta(mean * k, (1 - mean) * k, n)


def draw_inputs(rng, n, arm, p0_med, p_sub=None, emax_mult=1.0, sig_scale=1.0):
    """Sample all uncertain inputs for one arm."""
    a = ARMS[arm]
    s = sig_scale
    g = {}
    g["P0"] = _ln(rng, p0_med, PRIORS["P0"]["sig"] * s, n)
    g["k_deg"] = _ln(rng, PRIORS["k_deg"]["med"], PRIORS["k_deg"]["sig"] * s, n)
    g["S_post"] = _beta_ms(rng, PRIORS["S_post"]["mean"], PRIORS["S_post"]["sd"] * s, n)
    ps = PRIORS["p_sub"]["val"] if p_sub is None else p_sub
    sub = rng.random(n) < ps
    er0 = _ln(rng, PRIORS["ER0"]["med"], PRIORS["ER0"]["sig"] * s, n)
    g["ER0"] = np.where(sub, np.maximum(er0, 1.0), 1.0)
    g["A_m"] = _ln(rng, PRIORS["A_m"]["med"], PRIORS["A_m"]["sig"] * s, n)
    g["EC50"] = _ln(rng, PRIORS["EC50"]["med"], PRIORS["EC50"]["sig"] * s, n)
    g["Emax"] = _ln(rng, PRIORS["Emax"]["med"] * emax_mult, PRIORS["Emax"]["sig"] * s, n)
    g["h"] = _ln(rng, PRIORS["h"]["med"], PRIORS["h"]["sig"] * s, n)
    g["C_tox"] = _ln(rng, PRIORS["C_tox"]["med"], PRIORS["C_tox"]["sig"] * s, n)
    g["t_rev"] = _ln(rng, PRIORS["t_rev"]["med"], PRIORS["t_rev"]["sig"] * s, n)
    g["IC50"] = _ln(rng, PRIORS["IC50"]["med"], PRIORS["IC50"]["sig"] * s, n)
    # formulation-controlled (arm medians with dispersion)
    g["delta"] = _ln(rng, a["delta"], 0.30 * s, n)
    g["k_dil"] = _ln(rng, a["k_dil"], 0.30 * s, n)
    g["m"] = _ln(rng, 1.0, a["mism"] * s, n)
    g["Pr"] = _beta_ms(rng, max(a["prot"], 1e-3), 0.10 * s, n) if a["prot"] > 0 else np.zeros(n)
    g["k_dep"] = _ln(rng, a["k_dep"], 0.30 * s, n)
    g["k_rel"] = _ln(rng, a["k_rel"], 0.20 * s, n)
    return g


def integrate(g, arm, dt=DT, t_end=T_END, pe=None, inh=None, track=False):
    """Vectorised explicit-Euler integration.  Returns dict of outputs."""
    a = ARMS[arm]
    pe_dose = a["pe"] if pe is None else pe
    inh_dose = a["inh"] if inh is None else inh
    n = len(g["P0"])
    V = g["A_m"] * g["delta"]
    Pt = np.ones(n); P = np.zeros(n)
    Et = np.full(n, pe_dose); E = np.zeros(n)
    It = np.full(n, inh_dose); I = np.zeros(n)
    Abs = np.zeros(n); Deg = np.zeros(n); Dil = np.zeros(n)
    DH = np.zeros(n); t_last = np.zeros(n)
    kdeg_eff = g["k_deg"] * (1.0 - g["Pr"])
    steps = int(round(t_end / dt))
    traj = [] if track else None
    for i in range(steps):
        t = i * dt
        Ce = 1000.0 * E / V      # mmol/mL -> mM
        Ci = 1000.0 * I / V
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(Ce > 0, (Ce / g["EC50"]) ** g["h"], 0.0)
            H = x / (1.0 + x)
            y = np.where(Ce > 0, (Ce / g["C_tox"]) ** 3.0, 0.0)
            Hd = y / (1.0 + y)
        Peff = g["P0"] * (1.0 + g["Emax"] * H)
        ER = 1.0 + (g["ER0"] - 1.0) / (1.0 + Ci / g["IC50"])
        kabs = Peff / (g["delta"] * ER)
        rP = g["k_rel"] * Pt
        rE = g["k_rel"] * g["m"] * Et
        rI = g["k_rel"] * g["m"] * It
        a_flux = kabs * P
        d_flux = kdeg_eff * P
        w_flux = g["k_dil"] * P
        Pt = Pt - rP * dt
        P = P + (rP - a_flux - d_flux - w_flux) * dt
        Abs = Abs + a_flux * dt
        Deg = Deg + d_flux * dt
        Dil = Dil + w_flux * dt
        Et = Et - rE * dt
        E = E + (rE - (g["k_dil"] + g["k_dep"]) * E) * dt
        It = It - rI * dt
        I = I + (rI - g["k_dil"] * I) * dt
        DH = DH + Hd * dt
        t_last = np.where(Ce > 0.1 * g["EC50"], t + dt, t_last)
        if track:
            traj.append((t + dt, np.median(Ce / g["EC50"]), np.median(Ce / g["C_tox"]),
                         np.median(P), np.median(Abs)))
    F = Abs * g["S_post"]
    safe = (DH <= DH_MAX) & ((t_last + g["t_rev"]) <= T_REVERSE_MAX)
    out = dict(F=F, safe=safe, DH=DH, t_last=t_last, Abs=Abs, Deg=Deg, Dil=Dil,
               rem=Pt + P, closure=Pt + P + Abs + Deg + Dil)
    if track:
        out["traj"] = np.array(traj)
    return out


def run_arm(n, arm, p0_med, seed_offset=0, **kw):
    rng = np.random.default_rng(SEED + seed_offset)
    sub = {k: kw.pop(k) for k in ("pe", "inh", "dt") if k in kw}
    g = draw_inputs(rng, n, arm, p0_med, **kw)
    return g, integrate(g, arm, **sub)


def calibrate(n=8000):
    """Scale P0 median so that the reference arm R1 has median F = 1.2 %."""
    lo, hi = 1e-6, 1.0
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        _, o = run_arm(n, "R1_reference_PE", mid, seed_offset=1)
        if np.median(o["F"]) > F_BASELINE_MID:
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


def summarise(o):
    F = o["F"]; safe = o["safe"]
    ok = (F >= F_TARGET) & safe
    n = len(F)
    return dict(
        median=float(np.median(F)), p05=float(np.percentile(F, 5)), p95=float(np.percentile(F, 95)),
        p_F_ge_2=float(np.mean(F >= 0.02)), p_F_ge_4=float(np.mean(F >= F_TARGET)),
        p_F_ge_5=float(np.mean(F >= F_STRETCH)), p_safe=float(np.mean(safe)),
        p_success=float(np.mean(ok)),
        se_success=float(np.sqrt(np.mean(ok) * (1 - np.mean(ok)) / n)),
        median_DH=float(np.median(o["DH"])), median_t_clear=float(np.median(o["t_last"])),
    )


def spearman(x, y):
    """Spearman rank correlation with average ranks for ties (v2; v1 broke ties arbitrarily, which changed 18 reported
    sensitivity values by <= 0.03 between platforms)."""
    from scipy.stats import rankdata
    rx = rankdata(x).astype(float)
    ry = rankdata(y).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    return float((rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)))


def full_analysis(n=N_DEFAULT, outdir="results", figures=True):
    os.makedirs(outdir, exist_ok=True)
    res = dict(meta=dict(seed=SEED, n=n, dt=DT, t_end_h=T_END, F_target=F_TARGET,
                         DH_max_h=DH_MAX, reversal_max_h=T_REVERSE_MAX,
                         baseline_anchor=F_BASELINE_MID))
    p0 = calibrate()
    res["meta"]["P0_median_cm_per_h"] = p0
    res["priors"] = {k: {kk: vv for kk, vv in v.items()} for k, v in PRIORS.items()}
    res["priors"]["P0"]["med"] = p0
    res["arms_def"] = ARMS

    # 1) arms
    store = {}
    res["arms"] = {}
    for j, arm in enumerate(ARMS):
        g, o = run_arm(n, arm, p0, seed_offset=100 + j)
        store[arm] = (g, o)
        res["arms"][arm] = summarise(o)
        res["arms"][arm]["median_F_per_mmol_PE"] = (
            float(np.median(o["F"]) / ARMS[arm]["pe"]) if ARMS[arm]["pe"] > 0 else None)
        res["arms"][arm]["max_closure_error"] = float(np.max(np.abs(o["closure"] - 1.0)))

    # 2) efflux-substrate scenarios
    res["scenarios_p_sub"] = {}
    for ps in (0.0, 0.5, 1.0):
        res["scenarios_p_sub"][str(ps)] = {}
        for j, arm in enumerate(("R1_reference_PE", "R2_reference_PE_efflux", "P1_dual_lock", "P2_no_efflux_module")):
            _, o = run_arm(n // 2, arm, p0, seed_offset=200 + j, p_sub=ps)
            res["scenarios_p_sub"][str(ps)][arm] = summarise(o)

    # 3) prior-width robustness
    res["prior_width"] = {}
    for lab, sc in (("narrow_0.5x", 0.5), ("base_1.0x", 1.0), ("wide_1.5x", 1.5)):
        res["prior_width"][lab] = {}
        for j, arm in enumerate(("R1_reference_PE", "P1_dual_lock")):
            _, o = run_arm(n // 2, arm, p0, seed_offset=300 + j, sig_scale=sc)
            res["prior_width"][lab][arm] = summarise(o)

    # 4) global sensitivity (Spearman) for P1 on F, and on SUCCESS indicator
    g, o = store["P1_dual_lock"]
    ok = ((o["F"] >= F_TARGET) & o["safe"]).astype(float)
    sens = {}
    for k in ("P0", "k_deg", "S_post", "ER0", "A_m", "EC50", "Emax", "h", "C_tox", "t_rev",
              "IC50", "delta", "k_dil", "m", "Pr", "k_dep", "k_rel"):
        sens[k] = dict(rho_F=spearman(g[k], o["F"]), rho_success=spearman(g[k], ok))
    res["sensitivity_P1"] = sens

    # 5) design-window grids
    m = max(2000, n // 10)
    deltas = np.array([0.15, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00])
    doses = np.array([0.10, 0.20, 0.35, 0.50, 0.75, 1.00])
    G1 = np.zeros((len(doses), len(deltas)))
    for a_i, d in enumerate(doses):
        for b_i, dl in enumerate(deltas):
            ARMS["_grid"] = dict(ARMS["P1_dual_lock"]); ARMS["_grid"]["delta"] = float(dl)
            rng = np.random.default_rng(SEED + 400)
            gg = draw_inputs(rng, m, "_grid", p0)
            oo = integrate(gg, "_grid", pe=float(d))
            G1[a_i, b_i] = np.mean((oo["F"] >= F_TARGET) & oo["safe"])
    del ARMS["_grid"]
    res["grid_confinement_dose"] = dict(doses_mmol=doses.tolist(), delta_cm=deltas.tolist(), p_success=G1.tolist())

    emult = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    psubs = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    G2 = np.zeros((len(psubs), len(emult)))
    G2r = np.zeros((len(psubs), len(emult)))
    for a_i, ps in enumerate(psubs):
        for b_i, em in enumerate(emult):
            rng = np.random.default_rng(SEED + 500)
            gg = draw_inputs(rng, m, "P1_dual_lock", p0, p_sub=float(ps), emax_mult=float(em))
            oo = integrate(gg, "P1_dual_lock")
            G2[a_i, b_i] = np.mean((oo["F"] >= F_TARGET) & oo["safe"])
            rng = np.random.default_rng(SEED + 501)
            gg = draw_inputs(rng, m, "R1_reference_PE", p0, p_sub=float(ps), emax_mult=float(em))
            oo = integrate(gg, "R1_reference_PE")
            G2r[a_i, b_i] = np.mean((oo["F"] >= F_TARGET) & oo["safe"])
    res["grid_emax_psub"] = dict(emax_mult=emult.tolist(), p_sub=psubs.tolist(),
                                 p_success_P1=G2.tolist(), p_success_R1=G2r.tolist())

    # 6) go / no-go: which PE kinetic window is required (conditional on P1 architecture)
    #    probability of success among draws binned by Emax*EC50^-1 composite potency
    pot = g["Emax"] / g["EC50"]
    qs = np.quantile(pot, [0, .25, .5, .75, 1.0])
    gate = []
    for lo_, hi_ in zip(qs[:-1], qs[1:]):
        sel = (pot >= lo_) & (pot <= hi_)
        gate.append(dict(potency_lo=float(lo_), potency_hi=float(hi_),
                         p_success=float(np.mean(ok[sel])), n=int(sel.sum())))
    res["potency_quartiles_P1"] = gate

    # 6b) efflux-module dose-response (peptide always an efflux substrate; architecture P1)
    inh_doses = [0.0, 0.01, 0.033, 0.1, 0.3]
    res["inhibitor_dose_response"] = []
    for d in inh_doses:
        rng = np.random.default_rng(SEED + 700)
        gg = draw_inputs(rng, m, "P1_dual_lock", p0, p_sub=1.0)
        oo = integrate(gg, "P1_dual_lock", inh=d)
        res["inhibitor_dose_response"].append(dict(inh_mmol=d, **summarise(oo)))

    # 7) time course (median-parameter view of reversibility)
    tc = {}
    for arm in ("R1_reference_PE", "P1_dual_lock"):
        rng = np.random.default_rng(SEED + 600)
        gg = draw_inputs(rng, 4000, arm, p0)
        oo = integrate(gg, arm, track=True)
        tc[arm] = oo["traj"].tolist()
    res["timecourse"] = tc

    # determinism hash of the headline numbers
    head = json.dumps({k: res["arms"][k]["p_success"] for k in res["arms"]}, sort_keys=True)
    res["meta"]["headline_hash"] = hashlib.sha256(head.encode()).hexdigest()[:16]

    with open(os.path.join(outdir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    if figures:
        make_figures(res, store, outdir)
    return res


def make_figures(res, store, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 130, "savefig.dpi": 300})
    labels = {"R0_peptide_alone": "R0\npeptide\nalone", "R1_reference_PE": "R1\nreference\nPE",
              "R2_reference_PE_efflux": "R2\nref. PE +\nefflux", "P1_dual_lock": "P1\nDual-Lock",
              "P2_no_efflux_module": "P2\n− efflux", "P3_no_colocalisation": "P3\n− co-\nlocalisation",
              "P4_no_protection": "P4\n− protection", "P5_no_self_limiting_PE": "P5\n− self-\nlimiting PE"}
    arms = list(ARMS.keys())
    # Fig 2 : distributions
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    data = [np.clip(store[a][1]["F"] * 100, 1e-3, None) for a in arms]
    parts = ax.violinplot([np.log10(d) for d in data], showmedians=True, widths=0.8)
    for pc in parts["bodies"]:
        pc.set_alpha(0.55)
    ax.set_xticks(range(1, len(arms) + 1)); ax.set_xticklabels([labels[a] for a in arms], fontsize=7)
    yt = [0.01, 0.1, 1, 4, 10, 50]
    ax.set_yticks([np.log10(v) for v in yt]); ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.axhline(np.log10(4), ls="--", lw=1, color="C3"); ax.text(8.45, np.log10(4) + 0.05, "4 % target", color="C3", fontsize=7, ha="right")
    ax.axhline(np.log10(1.2), ls=":", lw=1, color="gray"); ax.text(8.45, np.log10(1.2) - 0.16, "1.2 % baseline anchor", color="gray", fontsize=7, ha="right")
    ax.set_ylabel("absolute oral bioavailability F (%)"); ax.set_ylim(-1.6, 2.0); ax.set_xlim(0.4, 8.6)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig2_distributions.png")); plt.close(fig)

    # Fig 3 : success probabilities & scenarios
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1.5, 1]})
    ps = [res["arms"][a]["p_success"] * 100 for a in arms]
    se = [res["arms"][a]["se_success"] * 100 * 1.96 for a in arms]
    axs[0].bar(range(len(arms)), ps, yerr=se, color=["#999", "#999", "#999", "#1f77b4"] + ["#9ecae1"] * 4, capsize=2)
    axs[0].set_xticks(range(len(arms))); axs[0].set_xticklabels([labels[a].replace("\n", " ") for a in arms], rotation=60, ha="right", fontsize=7)
    axs[0].set_ylabel("P(F ≥ 4 % and reversible), %")
    sc = res["scenarios_p_sub"]
    xs = np.arange(3); w = 0.2
    for i, arm in enumerate(("R1_reference_PE", "R2_reference_PE_efflux", "P1_dual_lock", "P2_no_efflux_module")):
        axs[1].bar(xs + (i - 1.5) * w, [sc[k][arm]["p_success"] * 100 for k in ("0.0", "0.5", "1.0")], w,
                   label=labels[arm].replace("\n", " "))
    axs[1].set_xticks(xs); axs[1].set_xticklabels(["never a\nsubstrate", "50 %", "always a\nsubstrate"], fontsize=7)
    axs[1].set_xlabel("peptide is an efflux substrate"); axs[1].set_ylim(0, 80); axs[1].legend(fontsize=5.5, frameon=False, loc="upper left", ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig3_success_scenarios.png")); plt.close(fig)

    # Fig 4 : tornado
    s = res["sensitivity_P1"]
    keys = sorted(s, key=lambda k: abs(s[k]["rho_F"]), reverse=True)[:10][::-1]
    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    ax.barh(keys, [s[k]["rho_F"] for k in keys], color=["#d62728" if s[k]["rho_F"] < 0 else "#1f77b4" for k in keys])
    ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("Spearman ρ (input vs F), architecture P1")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig4_tornado.png")); plt.close(fig)

    # Fig 5 : design windows
    g1 = res["grid_confinement_dose"]; g2 = res["grid_emax_psub"]
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.2))
    im = axs[0].imshow(np.array(g1["p_success"]), origin="lower", aspect="auto", vmin=0, vmax=1, cmap="viridis")
    axs[0].set_xticks(range(len(g1["delta_cm"]))); axs[0].set_xticklabels(g1["delta_cm"], fontsize=7)
    axs[0].set_yticks(range(len(g1["doses_mmol"]))); axs[0].set_yticklabels(g1["doses_mmol"], fontsize=7)
    axs[0].set_xlabel("niche thickness δ (cm)"); axs[0].set_ylabel("PE dose (mmol)")
    axs[0].set_title("P1: P(success)", fontsize=8)
    im2 = axs[1].imshow(np.array(g2["p_success_P1"]), origin="lower", aspect="auto", vmin=0, vmax=1, cmap="viridis")
    axs[1].set_xticks(range(len(g2["emax_mult"]))); axs[1].set_xticklabels(g2["emax_mult"], fontsize=7)
    axs[1].set_yticks(range(len(g2["p_sub"]))); axs[1].set_yticklabels(g2["p_sub"], fontsize=7)
    axs[1].set_xlabel("PE potency multiplier on E$_{max}$"); axs[1].set_ylabel("P(peptide is efflux substrate)")
    axs[1].set_title("P1: P(success)", fontsize=8)
    fig.colorbar(im2, ax=axs, fraction=0.025, pad=0.02)
    fig.savefig(os.path.join(outdir, "fig5_design_window.png"), bbox_inches="tight"); plt.close(fig)

    # Fig 6 : time course
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9))
    for arm, c in (("R1_reference_PE", "#777"), ("P1_dual_lock", "#1f77b4")):
        tr = np.array(res["timecourse"][arm])
        axs[0].plot(tr[:, 0], tr[:, 2], color=c, label=labels[arm].replace("\n", " "))
        axs[1].plot(tr[:, 0], tr[:, 4] * 100, color=c)
    axs[0].axhline(0.1, ls=":", color="k", lw=0.7)
    axs[0].set_xlabel("time (h)"); axs[0].set_ylabel("median PE exposure / C$_{tox}$"); axs[0].legend(fontsize=7, frameon=False)
    axs[1].set_xlabel("time (h)"); axs[1].set_ylabel("median cumulative apical uptake (% dose)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig6_timecourse.png")); plt.close(fig)


def self_tests():
    ok = True
    def chk(name, cond):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + name)
        ok = ok and bool(cond)
    p0 = calibrate(4000)
    g, o = run_arm(3000, "P1_dual_lock", p0, seed_offset=7)
    chk("mass balance closes to 1e-9 (peptide)", np.max(np.abs(o["closure"] - 1)) < 1e-9)
    chk("F within [0,1]", (o["F"] >= 0).all() and (o["F"] <= 1).all())
    g2, o2 = run_arm(3000, "P1_dual_lock", p0, seed_offset=7)
    chk("deterministic (same seed -> identical F)", np.array_equal(o["F"], o2["F"]))
    _, oa = run_arm(3000, "R1_reference_PE", p0, seed_offset=1)
    chk("calibration: R1 median F within 5% of 1.2 %", abs(np.median(oa["F"]) - F_BASELINE_MID) / F_BASELINE_MID < 0.05)
    _, ob = run_arm(3000, "R0_peptide_alone", p0, seed_offset=1)
    chk("PE raises median F over peptide alone", np.median(oa["F"]) > np.median(ob["F"]))
    gx = draw_inputs(np.random.default_rng(5), 2000, "P1_dual_lock", p0, p_sub=1.0)
    f1 = integrate(gx, "P1_dual_lock", inh=0.0)["F"]; f2 = integrate(gx, "P1_dual_lock", inh=0.033)["F"]
    chk("efflux module never lowers F when peptide is a substrate (same draws)", np.all(f2 >= f1 - 1e-12))
    gy = draw_inputs(np.random.default_rng(5), 2000, "P1_dual_lock", p0, p_sub=0.0)
    d1 = integrate(gy, "P1_dual_lock", inh=0.0)["F"]; d2 = integrate(gy, "P1_dual_lock", inh=0.033)["F"]
    chk("efflux module has no effect if peptide is not a substrate", np.allclose(d1, d2))
    gz = draw_inputs(np.random.default_rng(9), 1500, "P1_dual_lock", p0)
    fa = integrate(gz, "P1_dual_lock", dt=0.01)["F"]; fb = integrate(gz, "P1_dual_lock", dt=0.0025)["F"]
    chk("time-step convergence: median F changes < 2 % for dt/4", abs(np.median(fa) - np.median(fb)) / np.median(fb) < 0.02)
    chk("spearman of identical vectors = 1", abs(spearman(np.arange(50.0), np.arange(50.0)) - 1) < 1e-12)
    print("ALL PASS" if ok else "SOME FAILED")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=N_DEFAULT)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    if a.test:
        sys.exit(0 if self_tests() else 1)
    if a.run:
        r = full_analysis(a.n, a.out)
        print(json.dumps({k: round(v["p_success"], 3) for k, v in r["arms"].items()}, indent=1))
    if not (a.test or a.run):
        ap.print_help()
