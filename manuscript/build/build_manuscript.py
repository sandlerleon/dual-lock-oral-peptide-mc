# -*- coding: utf-8 -*-
"""Manuscript for the International Journal of Pharmaceutics, built from results_ext/*.json so that every number in the text is the model's.

    python build_manuscript.py <results_dir> <figures_dir> <out_dir>
"""
import json
import os
import re
import sys

import build_infra as BI
from refs_lib import LIB
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
import text_front

RES = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "../code/results_ext")
FIGDIR = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else os.path.join(RES, "figures"))
OUT = os.path.abspath(sys.argv[3] if len(sys.argv) > 3 else "../repo/manuscript")
HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(RES, "ext_%s.json" % n), encoding="utf-8"))
V1 = json.load(open(os.path.join(HERE, "results_v1.json"), encoding="utf-8"))
DOIS = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))
GITHUB = "https://github.com/sandlerleon/dual-lock-oral-peptide-mc"
TITLE = "Mechanistic Monte Carlo assessment of colocalized permeation enhancement and efflux modulation for oral peptide delivery"
EQ, PA, CV, CO, SO, PR, FA = J("equiv"), J("paired"), J("converge"), J("corr"), J("sobol"), J("pareto"), J("fals")
ST, AB, VA = J("struct"), J("abc"), J("valid")
AR = PA["arms"]


def p(x, d=1):
    return ("%." + str(d) + "f") % (100 * x)


def f(x, d=2):
    return ("%." + str(d) + "f") % x


def pp(x, d=1):
    return (("%+." + str(d) + "f") % (100 * x)).replace("-", "−")


def pct(x, d=1):
    return p(x, d) + "%"


B = list(text_front.B)
add = B.append
FIGN = {}
TABN = {2: "Table 2", 3: "Table 3"}
_fig = [0]


def fig(key, caption, w=6.3):
    _fig[0] += 1; FIGN[key] = "Fig. %d" % _fig[0]
    add(("fig", os.path.join(FIGDIR, key + ".png"), "**Fig. %d.** %s" % (_fig[0], caption), w))


def F(key):
    return FIGN[key]


def tbl(num, caption, header, rows, widths=None, size=15):
    TABN[num] = "Table %d" % num
    add(("table", "**Table %d.** %s" % (num, caption), header, rows, widths, size))


# Figure numbers are assigned in order of appearance; register the schematic first (it is placed in section 2 by insertion below).
FIGN.update({"fig1_schematic": "Fig. 1", "fig2_arms": "Fig. 2", "fig3_ablation": "Fig. 3", "fig4_sobol": "Fig. 4", "fig8_convergence_correlation": "Fig. 5", "fig5_pareto": "Fig. 6",
             "fig6_falsification": "Fig. 7", "fig11_mechanism_controls": "Fig. 8", "fig12_optimized_vs_nominal": "Fig. 9", "fig7_structural": "Fig. 10", "fig9_plausibility": "Fig. 11", "fig10_constrained": "Fig. 12"})
TABN_SEC = {"decomp": "Table 7", "optim": "Table 9"}
_fig[0] = 0


def fig_fixed(key, caption, w=6.3):
    _fig[0] += 1
    assert FIGN[key] == "Fig. %d" % _fig[0], (key, FIGN[key], _fig[0])
    add(("fig", os.path.join(FIGDIR, key + ".png"), "**%s.** %s" % (FIGN[key], caption), w))


# insert the schematic at the end of section 2 (after the "What the niche thickness is" paragraph)
_idx = max(i for i, b in enumerate(B) if b[0] == "p" and b[1].startswith("*What the niche thickness is.*"))
B.insert(_idx + 1, ("fig", os.path.join(FIGDIR, "fig1_schematic.png"), "**Fig. 1.** Structure of the Dual-Lock model. Tablet pools release peptide, enhancer and modulator into a co-localized niche of volume $V=A_m\\delta$; "
                    "uptake, proteolysis and wash-out compete for the peptide, and enhancer clearance and epithelial recovery set the duration of the exposure.", 6.2))
_fig[0] = 1

r1, p1 = AR["R1_reference_PE"], AR["P1_dual_lock"]
dom = PA["dominance_P1_vs_R1"]
vp1 = PA["paired_vs_P1"]; vr1 = PA["paired_vs_R1"]


def dsc(d):
    return "%s pp (95%% CI %s to %s)" % (pp(d[0]), pp(d[1]), pp(d[2]))


# --------------------------------------------------------------------------------------------------------------------- 5 results
add(("h1", "5. Results"))
add(("h2", "5.1 Reference comparison and sampler equivalence"))
zs = [EQ[k]["z_success"] for k in EQ]
add(("p", "With the reference arm calibrated to a median of 1.2%%, the Dual-Lock arm P1 has a median bioavailability of %s (5th to 95th percentile %s to %s) at half the enhancer dose, and %s of draws satisfy both the 4%% target and the reversibility rule, against %s for the reference (%s, %s). "
          "The enhancer-free peptide has a median of %s. The independent re-implementation of the sampler used for all paired analyses reproduced the original success probabilities within sampling error for four test arms (|z| $\\le$ %s)."
          % (pct(p1["median"]), pct(p1["p05"]), pct(p1["p95"]), pct(p1["p_success"]), pct(r1["p_success"]), TABN.get(4, "Table 4"), FIGN["fig2_arms"], pct(AR["R0_peptide_alone"]["median"]), f(max(abs(z) for z in zs), 1))))
arm_order = [("R0_peptide_alone", "R0 peptide alone"), ("R1_reference_PE", "R1 reference PE"), ("R2_reference_PE_efflux", "R2 R1 + efflux module"), ("P1_dual_lock", "P1 Dual-Lock"),
             ("P2_no_efflux_module", "P2 no efflux module"), ("P3_no_colocalisation", "P3 no co-localization"), ("P4_no_protection", "P4 reference protection"),
             ("P5_no_self_limiting_PE", "P5 no self-limiting PE"), ("M0_matrix_only", "M0 matrix only (negative control)")]
rows = []
for k, lab in arm_order:
    s = AR[k]; rows.append([lab, pct(s["median"], 2 if s["median"] < 0.02 else 1), "%s to %s" % (pct(s["p05"]), pct(s["p95"])), pct(s["p_F_ge_4"], 0), pct(s["p_safe"], 0), "%s ± %s" % (pct(s["p_success"], 0), p(1.96 * s["se"], 1))])
tbl(4, "Bioavailability and success by arm (paired draws, $n=20{,}000$ each). Success: $F\\ge4\\%$ and reversible. Probabilities are conditional on the priors and are not empirical probabilities.",
    ["Arm", "Median $F$", "5th to 95th percentile", "P($F\\ge4$%)", "P(reversible)", "P(success) $\\pm$ 95% MC error"], rows, [2.8, 1.0, 1.6, 1.0, 1.1, 1.8], 15)
fig_fixed("fig2_arms", "(a) Median and 5th to 95th percentile of absolute bioavailability by arm (log scale; dashed line, 4% target). (b) Probability of success with 95% Monte Carlo error. "
                       "The reference arm is calibrated to a median of 1.2%.", 6.4)

add(("h2", "5.2 Paired contrasts, ablations and negative controls"))
add(("p", "Because all arms share the same draws, differences can be read draw by draw. P1 gave a higher bioavailability than R1 in %s of draws and a lower damage-hours value in %s; it succeeded where R1 failed in %s of draws and failed where R1 succeeded in %s (median bioavailability ratio %s, median damage-hours ratio %s). "
          "The paired success difference was %s."
          % (pct(dom["p_F_higher"], 2), pct(dom["p_DH_lower"]), pct(dom["p_P1_succeeds_R1_fails"]), pct(dom["p_R1_succeeds_P1_fails"]), f(dom["median_ratio_F"], 1), f(dom["median_DH_ratio"], 2), dsc(vr1["P1_dual_lock"]["d_success"]))))
add(("p", "Removing one component from P1 (%s, left) shows where the benefit arises. Reverting the niche thickness to the reference value costs %s of success, reverting enhancer clearance (P5) %s, reverting the wash-out rate %s, removing the modulator (P2) %s, "
          "and reverting protection (P4) %s. Reverting the release mismatch changes success by only %s. Adding single components to R1 (right) reproduces the picture: the thin niche alone adds %s, but at the price of a lower probability of reversibility (%s against %s for R1); "
          "the whole confinement bundle without a self-limiting enhancer gives a median bioavailability of %s but only %s reversible draws, so success changes by %s, that is, it falls."
          % (FIGN["fig3_ablation"], dsc(vp1["P1_minus_delta"]["d_success"]), dsc(vp1["P5_no_self_limiting_PE"]["d_success"]), dsc(vp1["P1_minus_k_dil"]["d_success"]), dsc(vp1["P2_no_efflux_module"]["d_success"]), dsc(vp1["P4_no_protection"]["d_success"]),
             dsc(vp1["P1_minus_mism"]["d_success"]), dsc(vr1["R1_plus_delta"]["d_success"]), pct(AR["R1_plus_delta"]["p_safe"], 0), pct(r1["p_safe"], 0), pct(AR["R1_plus_confinement_bundle"]["median"]), pct(AR["R1_plus_confinement_bundle"]["p_safe"], 0),
             dsc(vr1["R1_plus_confinement_bundle"]["d_success"]))))
add(("p", "Two results differ from the original design intuition. First, the release rate assumed for the matrix (1.2 h$^{-1}$, slower than the reference 2.0 h$^{-1}$) lowers success: reverting it raises P(success) by %s points, so the model favors fast, synchronous release into the thin niche, "
          "not slow release. Second, the matrix-only negative control (P1 architecture, no enhancer and no modulator) has a median bioavailability of %s against %s for peptide alone, a gain of %s-fold from confinement alone. "
          "Published work found that absorption of semaglutide requires co-formulation with SNAC [@Buckley2018], so this gain is a consequence of the lumped structure (a smaller volume per contact area raises the peptide uptake rate $P_{eff}/\\delta$ and lowers dilution) and is not supported by data. "
          "It is flagged as a model assumption that the first experiment (Section 8) must test."
          % (f(100 * vp1["P1_minus_k_rel"]["d_success"][0], 1), pct(AR["M0_matrix_only"]["median"], 2), pct(AR["R0_peptide_alone"]["median"], 2), f(AR["M0_matrix_only"]["median"] / AR["R0_peptide_alone"]["median"], 0))))
fig_fixed("fig3_ablation", "Paired changes in P(success) (percentage points, 95% CI) when one component is removed from Dual-Lock P1 (left) or added to the reference R1 (right). "
                           "All arms use the same random draws.", 6.4)

abl = [("P1_minus_delta", "P1 with reference niche thickness", "paired_vs_P1"), ("P1_minus_k_dil", "P1 with reference wash-out", "paired_vs_P1"), ("P1_minus_mism", "P1 with reference release mismatch", "paired_vs_P1"),
       ("P1_minus_k_rel", "P1 with reference release rate (2.0 h$^{-1}$)", "paired_vs_P1"), ("P4_no_protection", "P1 with reference protection (P4)", "paired_vs_P1"), ("P5_no_self_limiting_PE", "P1 without self-limiting enhancer (P5)", "paired_vs_P1"),
       ("P2_no_efflux_module", "P1 without modulator (P2)", "paired_vs_P1"), ("P3_no_colocalisation", "P1 without co-localization (P3)", "paired_vs_P1"), ("P1_full_dose", "P1 at full enhancer dose (1.0 mmol)", "paired_vs_P1"),
       ("M0_matrix_only", "M0 matrix only, no enhancer or modulator", "paired_vs_P1"),
       ("R1_plus_delta", "R1 with thin niche (0.30 cm)", "paired_vs_R1"), ("R1_plus_k_dil", "R1 with slow wash-out (0.8 h$^{-1}$)", "paired_vs_R1"), ("R1_plus_k_rel", "R1 with slower release (1.2 h$^{-1}$)", "paired_vs_R1"),
       ("R1_plus_prot", "R1 with protection 0.70", "paired_vs_R1"), ("R1_plus_k_dep", "R1 with fast enhancer clearance (3 h$^{-1}$)", "paired_vs_R1"), ("R1_plus_confinement_bundle", "R1 with confinement bundle", "paired_vs_R1"),
       ("R2_reference_PE_efflux", "R1 with efflux module (R2)", "paired_vs_R1"), ("R1_half_dose", "R1 at half enhancer dose (0.5 mmol)", "paired_vs_R1")]
rows = []
for k, lab, rel in abl:
    s_ = AR[k]; d = PA[rel][k]["d_success"]
    rows.append([lab, pct(s_["median"], 2 if s_["median"] < 0.02 else 1), pct(s_["p_safe"], 0), pct(s_["p_success"], 0), "%s (%s to %s)" % (pp(d[0]), pp(d[1]), pp(d[2])), "P1" if rel == "paired_vs_P1" else "R1"])
tbl(5, "Ablations and negative controls on shared draws. $\Delta$ is the paired change in P(success) in percentage points (95% CI) relative to the arm named in the last column.",
    ["Modification", "Median $F$", "P(reversible)", "P(success)", "$\Delta$ (95% CI)", "vs"], rows, [3.6, 0.9, 1.0, 1.0, 2.0, 0.5], 14)

add(("h2", "5.3 Numerical and Monte Carlo convergence"))
rp1, rr1 = CV["replicates"]["P1_dual_lock"], CV["replicates"]["R1_reference_PE"]
dts = CV["dt"]; hz = CV["horizon"]["P1_dual_lock"]
add(("p", "Twenty independent replicate runs of 20,000 draws gave a P(success) for P1 of %s (sd %s percentage points, range %s to %s) and a median bioavailability of %s (sd %s points); for R1 the values were %s (sd %s points). "
          "Running estimates were stable from about $10^4$ draws and the batch-means standard error at $10^5$ draws was %s points (%s). Reducing the time step from 0.04 h to 0.0025 h changed the median bioavailability of P1 by %s, and from 0.01 h (the step used) to 0.0025 h by %s; "
          "success changed from %s to %s. Extending the horizon from 8 h to 24 h changed the median by %s (Supplementary Fig. S1)."
          % (pct(rp1["p_success_mean"]), f(100 * rp1["p_success_sd"], 2), p(rp1["p_success_range"][0]), p(rp1["p_success_range"][1]), pct(rp1["median_F_mean"], 2), f(100 * rp1["median_F_sd"], 2), pct(rr1["p_success_mean"]), f(100 * rr1["p_success_sd"], 2),
             f(100 * CV["batch_means_P1_dual_lock"]["se"], 2), FIGN["fig8_convergence_correlation"], p(abs(dts["0.04"]["median_F"] / dts["0.0025"]["median_F"] - 1)) + "%", p(abs(dts["0.01"]["median_F"] / dts["0.0025"]["median_F"] - 1)) + "%",
             pct(dts["0.01"]["p_success"]), pct(dts["0.0025"]["p_success"]), p(abs(hz["24.0"]["median_F"] / hz["8.0"]["median_F"] - 1)) + "%")))

add(("h2", "5.4 Global sensitivity and uncertainty decomposition"))
sp = SO["P1_dual_lock"]["log10_F"]["params"]; sg = SO["P1_dual_lock"]["log10_F"]["groups"]; ss = SO["P1_dual_lock"]["success"]["params"]; sgs = SO["P1_dual_lock"]["success"]["groups"]
top = sorted(sp, key=lambda k: -sp[k]["ST"])[:5]
NM = {"delta": "the niche thickness $\\delta$", "P0": "the baseline permeability $P_0$", "EC50": "the enhancer $EC_{50}$", "A_m": "the contact area $A_m$", "h": "the Hill coefficient", "Emax": "$E_{max}$", "S_post": "$S_{post}$", "k_deg": "$k_{deg}$",
      "k_dep": "the enhancer clearance $k_{dep}$", "t_rev": "the recovery time $t_{rev}$", "C_tox": "$C_{tox}$", "k_dil": "the wash-out rate"}
add(("p", "The variance of $\\log_{10}F$ for P1 is almost additive (sum of first-order indices %s). The largest total-order indices are %s. Grouped by origin, the total-order indices are %s for peptide and epithelium, %s for enhancer pharmacology, %s for formulation architecture and %s for the efflux modulator (%s). "
          "The uncertainty of the predicted bioavailability is thus shared about equally among peptide biology, enhancer pharmacology and the architecture, and the modulator potency contributes almost none, although removing the modulator costs %s of success on average (Section 5.2): "
          "the modulator shifts the mean for substrate peptides without adding variance of its own (the first-order index of the substrate indicator is %s)."
          % (f(sum(sp[k]["S1"] for k in sp), 2), ", ".join("%s (%s)" % (NM.get(k, k), f(sp[k]["ST"], 2)) for k in top), f(sg["peptide and epithelium"]["total"], 2), f(sg["enhancer pharmacology"]["total"], 2), f(sg["formulation architecture"]["total"], 2),
             f(sg["efflux modulator"]["total"], 2), FIGN["fig4_sobol"], dsc(vp1["P2_no_efflux_module"]["d_success"]), f(sp["sub"]["S1"], 2))))
tops = sorted(ss, key=lambda k: -ss[k]["ST"])[:6]
add(("p", "The success indicator is strongly interacting (sum of first-order indices %s). Its largest total-order indices are %s; group totals are %s (peptide and epithelium), %s (enhancer pharmacology), %s (architecture) and %s (modulator). "
          "The design lever with the largest effect on the *mean* (the niche, Section 5.2) therefore does not dominate the *uncertainty*: delivering the architecture is necessary, but the outcome still depends on the peptide and the enhancer."
          % (f(sum(ss[k]["S1"] for k in ss), 2), ", ".join("%s (%s)" % (NM.get(k, k), f(ss[k]["ST"], 2)) for k in tops), f(sgs["peptide and epithelium"]["total"], 2), f(sgs["enhancer pharmacology"]["total"], 2), f(sgs["formulation architecture"]["total"], 2), f(sgs["efflux modulator"]["total"], 2))))
fig_fixed("fig4_sobol", "Global variance-based sensitivity of Dual-Lock P1 ($N=2^{13}$ base samples, Saltelli and Jansen estimators, 95% bootstrap intervals). (a) First-order and total-order indices for $\\log_{10}F$ (twelve largest). "
                        "(b) Total-order indices of the four input groups for $\\log_{10}F$ and for the success indicator, P1 and R1.", 6.4)

add(("h2", "5.5 Correlated priors"))
sc = CO["scenarios"]; sw = CO["sweep_rho_EC50_Ctox"]; so = CO["sweep_other"]
allv = [v[a]["p_success"] for v in list(sc.values()) + list(sw.values()) + list(so.values()) for a in ("R1_reference_PE", "P1_dual_lock")]
adv = [k for k in so if k.startswith("adverse")][0]; fav = [k for k in so if k.startswith("favourable")][0]
add(("p", "With the assumed correlation set, P(success) for P1 changed from %s to %s and for R1 from %s to %s. Varying the coupling between enhancer potency and toxicity from 0 to 0.9 changed P(success) by at most %s percentage points, because the reversibility rule is governed by the recovery time and the enhancer clearance, not by damage-hours, "
          "in R1 and P1 (median damage-hours %s h and %s h against the 0.5 h limit). The adverse bundle of correlations gave %s (P1) and %s (R1), and the favorable bundle %s and %s (%s). The conclusion that P1 exceeds R1 does not depend on the assumed couplings."
          % (pct(sc["independent"]["P1_dual_lock"]["p_success"]), pct(sc["base_correlated"]["P1_dual_lock"]["p_success"]), pct(sc["independent"]["R1_reference_PE"]["p_success"]), pct(sc["base_correlated"]["R1_reference_PE"]["p_success"]),
             f(100 * (max(sw[r]["P1_dual_lock"]["p_success"] for r in sw) - min(sw[r]["P1_dual_lock"]["p_success"] for r in sw)), 2), f(r1["median_DH"], 3), f(p1["median_DH"], 3),
             pct(so[adv]["P1_dual_lock"]["p_success"]), pct(so[adv]["R1_reference_PE"]["p_success"]), pct(so[fav]["P1_dual_lock"]["p_success"]), pct(so[fav]["R1_reference_PE"]["p_success"]), FIGN["fig8_convergence_correlation"])))
fig_fixed("fig8_convergence_correlation", "Monte Carlo convergence and prior coupling. (a) Running P(success) of R1 and P1 with 95% error bands. (b) Running median bioavailability. (c) P(success) against the assumed rank coupling of enhancer potency and toxicity, $\\rho(EC_{50},C_{tox})$.", 6.4)

add(("h2", "5.6 Benefit-risk (Pareto) analysis"))
fam = PR["families"]; fR = fam["R1 reference (SNAC-like)"]; fP = fam["P1 Dual-Lock"]; f5 = fam["P5 Dual-Lock without self-limiting PE"]; f3 = fam["P3 Dual-Lock without co-localisation"]
at = lambda rows, d: [s for s in rows if abs(s["dose"] - d) < 1e-9][0]
best = max(fP, key=lambda s: s["p_success"]); bestR = max(fR, key=lambda s: s["p_success"]); best5 = max(f5, key=lambda s: s["p_success"]); best3 = max(f3, key=lambda s: s["p_success"])
add(("p", "A genuine benefit-risk improvement requires more than a higher bioavailability, so the enhancer dose was swept for the reference, Dual-Lock and two ablations and P($F\\ge4$%%) was plotted against the probability that the exposure is not reversible (%s). "
          "At 0.5 mmol P1 reaches P($F\\ge4$%%) of %s with %s of draws non-reversible; the reference needs 3 mmol to reach %s and then %s of its draws are non-reversible. At a similar risk (%s for R1 at 1.0 mmol, %s for P1 at 0.5 mmol) the reference reaches %s against %s for P1. "
          "The curve of P1 lies above that of R1 at every dose tested, so Dual-Lock dominates the reference in this plane. P(success) peaks at %s mmol for P1 (%s; %s at the 0.5 mmol design dose), saturates at %s for the reference (maximum at %g mmol), and peaks at only %s for P5 (%g mmol), "
          "which shows directly that a persistent enhancer raises bioavailability at the cost of reversibility. Removing co-localization (P3) leaves the exposure low but needs 3 mmol to reach %s."
          % (FIGN["fig5_pareto"], pct(at(fP, 0.5)["p_F_ge_4"], 0), pct(1 - at(fP, 0.5)["p_safe"], 0), pct(at(fR, 3.0)["p_F_ge_4"], 0), pct(1 - at(fR, 3.0)["p_safe"], 0), pct(1 - at(fR, 1.0)["p_safe"], 0), pct(1 - at(fP, 0.5)["p_safe"], 0), pct(at(fR, 1.0)["p_F_ge_4"], 0), pct(at(fP, 0.5)["p_F_ge_4"], 0),
             "%g" % best["dose"], pct(best["p_success"], 0), pct(at(fP, 0.5)["p_success"], 0), pct(bestR["p_success"], 0), bestR["dose"], pct(best5["p_success"], 0), best5["dose"], pct(best3["p_success"], 0))))
fig_fixed("fig5_pareto", "Benefit-risk analysis. (a) P($F\\ge4$%) against the probability that the exposure is not reversible for enhancer doses from 0.05 to 3 mmol (annotated: 0.5 and 1.0 mmol). "
                         "(b) Individual draws: bioavailability against damage-hours (dashed lines: 4% and 0.5 h).", 6.4)

add(("h2", "5.7 Falsification regions"))
d_list, dose_list = FA["deltas"], FA["doses"]
Gs = FA["grid_delta_dose"]["p_success"]; Gm = FA["grid_delta_dose"]["median"]; Gsafe = FA["grid_delta_dose"]["p_safe"]
j03 = d_list.index(0.30); j06 = d_list.index(0.60); j10 = d_list.index(1.00); j015 = d_list.index(0.15)
mx = FA["max_p_success_by_delta"]
add(("p", "The predefined benefit region B is the set of niche thickness and dose combinations for which the median bioavailability is at least 4%% and at least 80%% of draws are reversible (%s). Region B is non-empty for niche thicknesses up to %g cm. "
          "Along $\\delta$, the best P(success) over the dose range falls from %s at 0.15 cm through %s at 0.30 cm and %s at 0.60 cm to %s at 1.0 cm, with a best dose of %g mmol at 0.30 cm; the largest thickness at which some dose still achieves P(success) $\\ge$ 30%% is %g cm. "
          "Dual-Lock therefore fails the predefined benefit target if the matrix cannot confine the effective niche to %g cm or less (a %s-fold reduction of the reference volume per contact area), and it fails even the weaker 30%% target above %g cm (%s-fold). "
          "The second map shows that slowing the wash-out helps less than thinning the niche, and that in very thin niches a slow wash-out lowers success through persistent exposure (%s of draws at 0.05 cm with a wash-out of 0.2 h$^{-1}$ against %s at 4 h$^{-1}$)."
          % (FIGN["fig6_falsification"], FA["delta_max_for_B"], pct(mx[j015], 0), pct(mx[j03], 0), pct(mx[j06], 0), pct(mx[j10], 0), FA["best_dose_by_delta"][j03], FA["delta_max_for_p30"], FA["delta_max_for_B"], f(1.0 / FA["delta_max_for_B"], 1), FA["delta_max_for_p30"], f(1.0 / FA["delta_max_for_p30"], 1),
             pct(FA["grid_delta_kdil"]["p_success"][0][0], 0), pct(FA["grid_delta_kdil"]["p_success"][-1][0], 0))))
fig_fixed("fig6_falsification", "Falsification maps for Dual-Lock P1. (a) P(success), %, against niche thickness and enhancer dose; the white line encloses the benefit region B (median $F\\ge4$% and P(reversible) $\\ge0.80$). "
                                "(b) P(success) against niche thickness and wash-out rate at 0.5 mmol.", 6.4)

rows = []
for j, dl in enumerate(d_list):
    bi = max(range(len(dose_list)), key=lambda i: Gs[i][j])
    rows.append(["%g" % dl, "%g" % dose_list[bi], pct(Gs[bi][j], 0), pct(Gm[bi][j]), pct(Gsafe[bi][j], 0), "yes" if any(FA["region_B"][i][j] for i in range(len(dose_list))) else "no"])
tbl(6, "Falsification map in numbers: for each niche thickness $\delta$ (cm), the dose (mmol) with the highest P(success) and the corresponding median $F$ and P(reversible); the last column states whether benefit region B is non-empty. A best dose at the lower or upper edge of the grid (0.05 or 2 mmol) means that doses beyond it were not tested.",
    ["$\delta$ (cm)", "Best dose (mmol)", "Best P(success)", "Median $F$ at best dose", "P(reversible) at best dose", "Region B non-empty"], rows, [1.0, 1.2, 1.2, 1.6, 1.6, 1.4], 15)

# --------------------------------------------------------------------------------------------------------------------- 6 robustness
exec(open(os.path.join(HERE, "new_sections5.py"), encoding="utf-8").read())
add(("h1", "6. Robustness and falsification"))
add(("h2", "6.1 Structural uncertainty: stagnant layer and spatial transport"))
V = ST["variants"]
lum = V["lumped (base)"]; lay = V["lumped + stagnant layer"]; layN = V["lumped + stagnant layer, P0 not recalibrated"]
dtimes = ST["diffusion_times_h"]
add(("p", "The lumped niche assumes instantaneous mixing. With the diffusivity of a 3-kDa peptide in a gel of $10^{-6}$ cm$^2$ s$^{-1}$, the time for diffusion across 0.30 cm is %s h and across 1.0 cm %s h, longer than the 8 h horizon, so an unmixed niche would behave very differently. "
          "Table 10 and %s compare the lumped model with a serial stagnant layer (thickness 20 to 830 $\\mu$m, median 150 $\\mu$m [@Atuma2001], diffusivity in mucus [@Cu2009]) and with the spatial model for effective dispersion multipliers from 3 to $10^6$ times the molecular diffusivity. "
          "Each variant is recalibrated so that R1 again has a median of 1.2%%."
          % (f(dtimes["P1_peptide"], 0), f(dtimes["R1_peptide"], 0), FIGN["fig7_structural"])))
rows = []
LBL = lambda k: k.replace("spatial, dispersion x1e+06", "spatial, well-mixed limit (dispersion x10^6)").replace("spatial, dispersion x", "spatial, dispersion ×")
for k, v in V.items():
    if "not recalibrated" in k:
        continue
    if v.get("calibratable", True) is False:
        rows.append([LBL(k), "n/a", "n/a", "n/a", "n/a", "R1 cannot reach 1.2%% (maximum %s%% at $P_0$ = 40 cm h$^{-1}$)" % ("%.2g" % (100 * v["max_median_F_R1_at_P0_40"]))]); continue
    rows.append([LBL(k), "%sx" % f(v["p0_over_base"], 2), pct(v["R1_reference_PE"]["median"], 2), pct(v["P1_dual_lock"]["median"]), "%s / %s" % (pct(v["R1_reference_PE"]["p_success"], 0), pct(v["P1_dual_lock"]["p_success"], 0)),
                 "P3 %s, P5 %s" % (pct(v["P3_no_colocalisation"]["p_success"], 0), pct(v["P5_no_self_limiting_PE"]["p_success"], 0))])
tbl(10, "Structural uncertainty. Each variant is recalibrated to a median $F$ of 1.2% for R1 ($P_0$ relative to the lumped value is shown; a value far above 1 is implausible for a hydrophilic peptide). "
       "All variants use the same 2,000 draws per arm.",
    ["Structural model", "$P_0$/base", "R1 median $F$", "P1 median $F$", "P(success) R1 / P1", "Ablations: P(success)"], rows, [2.4, 0.8, 1.0, 1.0, 1.4, 2.4], 14)
fig_fixed("fig7_structural", "Structural uncertainty. (a) P(success) of R1 and P1 under the lumped model, the lumped model with a stagnant layer and the spatial diffusion-reaction model at several effective dispersion multipliers; annotations give the recalibrated $P_0$ relative to the base value. "
                             "(b) Median bioavailability of both arms (dashed line, 4%).", 6.4)

exec(open(os.path.join(HERE, "build_manuscript2.py"), encoding="utf-8").read())
os.makedirs(OUT, exist_ok=True)
exec(open(os.path.join(HERE, "build_manuscript3.py"), encoding="utf-8").read())

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    cit = BI.Harvard(LIB)
    BI.scan_blocks(cit, B)
    unused = [k for k in LIB if k not in cit.used]
    print("unused refs:", unused)
    cit.lib = {k: v for k, v in LIB.items() if k in cit.used}
    out = os.path.join(OUT, "DualLock_IJP_Manuscript.docx")
    BI.build(out, B, cit)
    print("wrote", out)
