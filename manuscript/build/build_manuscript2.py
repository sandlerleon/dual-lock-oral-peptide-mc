# -*- coding: utf-8 -*-
"""Part 2 of the manuscript builder (sections 6.1 interpretation to the declarations). Executed by build_manuscript.py in its namespace."""
import math
from scipy.stats import norm as _norm

vk = [k for k in V if "not recalibrated" not in k and V[k].get("calibratable", True) is not False and k.startswith("spatial")]
nc = [k for k in V if V[k].get("calibratable", True) is False]
v100 = V["spatial, dispersion x100"]; v1000 = V["spatial, dispersion x1000"]; v6 = V["spatial, dispersion x1e+06"]
nc_lbl = ", ".join(k.replace("spatial, dispersion ", "").replace("x", "×") for k in nc)


def tail(ratio):
    return 100 * (1 - _norm.cdf(math.log(ratio) / 0.6))


add(("p", "The structural variants lead to four findings. (i) With a stagnant layer in series and the baseline recalibrated, the median bioavailability of P1 is %s (lumped %s) and P(success) %s against %s for R1 (lumped %s and %s): the layer caps the enhanced permeability of the reference more than that of the confined arm. "
          "The recalibrated baseline permeability is %s times the lumped value, a value that the prior ($\\sigma_{\\ln}=0.6$) assigns an upper-tail probability of %s%%. "
          "(ii) In the spatial model the reference arm cannot reach the 1.2%% anchor when the effective dispersion is %s times the molecular diffusivity (maximum median %s, %s and %s). A conventional tablet in an unmixed 1 cm niche would absorb far less than observed, "
          "so the published baseline itself implies substantial mixing, and the lumped assumption is the one the anchor supports. "
          "(iii) Where the anchor can be met the Dual-Lock advantage does not disappear. At a dispersion multiplier of 100 the recalibrated $P_0$ is %s times the lumped value (tail probability %s%%, so this structure is implausible), P1 reaches a median of %s and P(success) of %s against %s for R1; "
          "at 1000 ($P_0$ %s times) the values are %s and %s against %s, and in the well-mixed limit they return to the lumped result (%s against %s). "
          "(iv) The lumped model (which equals the well-mixed limit of the spatial model) gives the smallest advantage of all calibratable structures, because transport resistance penalizes the thick-niche reference more than the thin-niche architecture; the lumped result is therefore not optimistic with respect to these alternatives. "
          "The high-advantage results of the poorly mixed structures should not be read as predictions, because they depend on an inflated baseline permeability and on a reference arm that is limited by diffusion rather than by permeability. "
          "Keeping $P_0$ at the lumped value instead of recalibrating gives reference bioavailabilities of %s (stagnant layer) and %s (dispersion $\\times100$), far below the anchor, which is a second way of seeing that the data exclude an unmixed niche."
          % (pct(lay["P1_dual_lock"]["median"]), pct(lum["P1_dual_lock"]["median"]), pct(lay["P1_dual_lock"]["p_success"], 0), pct(lay["R1_reference_PE"]["p_success"], 0), pct(lum["P1_dual_lock"]["p_success"], 0), pct(lum["R1_reference_PE"]["p_success"], 0),
             f(lay["p0_over_base"], 1), f(tail(lay["p0_over_base"]), 1), nc_lbl, "%.2g%%" % (100 * V[nc[0]]["max_median_F_R1_at_P0_40"]), "%.2g%%" % (100 * V[nc[1]]["max_median_F_R1_at_P0_40"]), "%.2g%%" % (100 * V[nc[2]]["max_median_F_R1_at_P0_40"]),
             f(v100["p0_over_base"], 1), f(tail(v100["p0_over_base"]), 1), pct(v100["P1_dual_lock"]["median"]), pct(v100["P1_dual_lock"]["p_success"], 0), pct(v100["R1_reference_PE"]["p_success"], 0),
             f(v1000["p0_over_base"], 1), pct(v1000["P1_dual_lock"]["median"]), pct(v1000["P1_dual_lock"]["p_success"], 0), pct(v1000["R1_reference_PE"]["p_success"], 0), pct(v6["P1_dual_lock"]["p_success"], 0), pct(v6["R1_reference_PE"]["p_success"], 0),
             pct(layN["R1_reference_PE"]["median"], 2), pct(V["spatial, dispersion x100, P0 not recalibrated"]["R1_reference_PE"]["median"], 2))))
add(("p", "The qualitative ordering is the same in every calibratable structure: P1 above R1, and P3 (no co-localization) and P5 (no self-limiting enhancer) well below P1. The quantitative advantage is structure-dependent. "
          "Whether the matrix in vivo provides the mixing assumed here, and how thick its stagnant layer is, are therefore measurable quantities on which the prediction depends (Section 8)."))

add(("h2", "6.2 Observation-conditioned analysis and comparison with published data"))
ak = [k for k in AB if k != "n"]
tbl(11, "Observation-conditioned analysis (rejection ABC, $n=10^5$ paired draws): outcome of Dual-Lock P1 among the draws in which the reference arm R1 falls in the stated published range.",
    ["Condition on R1", "Accepted draws", "P1 median $F$", "P1 5th to 95th percentile", "P1 P($F\\ge4$%)", "P1 P(success)"],
    [[k.replace("prior (no conditioning)", "none (prior)").replace(" %", "%"), "{:,}".format(AB[k]["n_accepted"]), pct(AB[k]["P1"]["median"]), "%s to %s" % (pct(AB[k]["P1"]["p05"]), pct(AB[k]["P1"]["p95"])), pct(AB[k]["P1"]["p_F_ge_4"], 0), pct(AB[k]["P1"]["p_success"], 0)] for k in ak],
    [2.8, 1.1, 1.0, 1.9, 1.0, 1.0], 14)
A2 = AB["R1 F in [0.4, 2.0] %"]; A1 = AB["R1 F in [0.4, 1.0] % (oral-semaglutide-like)"]; A3 = AB["R1 F in [0.8, 1.6] %"]
add(("p", "Table 11 shows that the Dual-Lock result depends on how well the reference performs. If R1 lies anywhere in 0.4 to 2.0%% (%s of the draws), P1 has a median of %s and P(success) of %s, close to the unconditioned value. "
          "If the reference performs like oral semaglutide, 0.4 to 1.0%% (%s of the draws), the P1 median falls to %s and P(success) to %s, about half the headline value, because the 4%% target is then much harder to reach from a lower starting point. "
          "Restricting R1 to 0.8 to 1.6%% around the anchor gives %s and %s. The headline probability therefore rests on the assumption that the reference tablet performs at the 1.2%% mid-point of the 1 to 2%% range; "
          "the conditioning shifts the baseline-determining parameters only slightly (medians of $P_0$ and $E_{max}$ by factors of %s and %s for the 0.4 to 1.0%% range) but not the ordering of the arms."
          % (pct(A2["frac"], 0), pct(A2["P1"]["median"]), pct(A2["P1"]["p_success"], 0), pct(A1["frac"], 0), pct(A1["P1"]["median"]), pct(A1["P1"]["p_success"], 0), pct(A3["P1"]["median"]), pct(A3["P1"]["p_success"], 0), f(A1["shift"]["P0"], 2), f(A1["shift"]["Emax"], 2))))
rv = VA["R1"]; pv = VA["P1"]; q = VA["qualitative"]; ob = VA["observed_equiv_sd_ln"]
add(("p", "Four comparisons with published data that were not used to calibrate the model give a mixed picture (%s). "
          "(a) *Spread.* The within-subject variability of oral semaglutide bioavailability is 137%%, equivalent to a standard deviation of the natural logarithm of %s [@Overgaard2021]; the reference arm has %s, because it combines variability with epistemic uncertainty about peptides, so the spread is wider than that of one product, as expected. "
          "(b) *Upper tail.* %s of R1 draws and %s of P1 draws exceed 10%%, %s of P1 draws exceed 20%%, and %s of P1 draws exceed 6%%, whereas the highest bioavailability reported for an enhancer-based oral peptide in the sources used is about 6%% in dogs for an engineered peptide [@Pechenov2021] and typical values are about 1%% [@Brayden2020]. "
          "The upper tail of P1 (95th percentile %s) lies beyond any published value; it is the upper end of the prior, not a prediction. "
          "(c) *Qualitative dependences.* The model predicts that bioavailability falls with a larger niche volume and with faster wash-out (median R1 bioavailability %s, %s and %s for $\\delta$ of 0.5, 1 and 2 cm; %s, %s and %s for wash-out multiplied by 0.5, 1 and 2), as oral semaglutide bioavailability falls with a larger water volume [@Overgaard2021]. "
          "(d) *Enhancer dose.* The model's bioavailability of R1 rises monotonically with the enhancer dose (%s at 0.5 mmol, %s at 2 mmol and %s at 8 mmol), whereas exposure to oral semaglutide was highest with 300 mg SNAC among 150 to 600 mg [@Granhall2019]. "
          "The model has no mechanism for such an optimum (for example a limit on local concentration, viscosity or pH); this structural deficiency matters for any dose-optimization claim. The transience reported by [!@Buckley2018] is reproduced qualitatively: "
          "%s of the final uptake of R1 occurs within 1 h and %s within 2 h."
          % (FIGN["fig9_plausibility"], f(ob, 2), f(rv["sd_lnF"], 2), pct(rv["frac_gt_10"], 1), pct(pv["frac_gt_10"], 0), pct(pv["frac_gt_20"], 0), pct(AB["prior (no conditioning)"]["P1"]["p_F_ge_4"] * 0 + pv["frac_gt_6"], 0), pct(p1["p95"], 0),
             pct(q["R1 delta=0.5 cm"], 2), pct(q["R1 delta=1 cm"], 2), pct(q["R1 delta=2 cm"], 2), pct(q["R1 k_dil x0.5"], 2), pct(q["R1 k_dil x1"], 2), pct(q["R1 k_dil x2"], 2),
             pct(q["R1 PE dose 0.5 mmol"], 2), pct(q["R1 PE dose 2 mmol"], 2), pct(q["R1 PE dose 8 mmol"], 1), pct(VA["fraction_absorbed_by_t"]["R1_reference_PE"]["1.0"], 0), pct(VA["fraction_absorbed_by_t"]["R1_reference_PE"]["2.0"], 0))))
fig_fixed("fig9_plausibility", "Calibration-independent checks. (a) Distribution of bioavailability for R1 and P1 (first 3,000 paired draws) with published reference values: oral semaglutide (0.4 to 1.0%, line at 0.8%), an engineered peptide in dogs (about 6%) and the 4% target. "
                               "(b) P(success) of P1 among the draws accepted under each conditioning rule.", 6.4)

exec(open(os.path.join(HERE, "new_sections6.py"), encoding="utf-8").read())
add(("h2", "6.3 Predefined failure conditions"))
nm_ = FA["delta_max_for_p30"]; nb_ = FA["delta_max_for_B"]
fail = [["Co-localization (Gate G1)", "effective niche thickness $\\delta>$ %g cm (less than 1.7-fold confinement of the reference 1.0 cm)" % nm_, "best P(success) over doses 0.05 to 2 mmol below 30% (Section 5.7)"],
        ["Benefit region B", "median $F<4$%% or P(reversible) $<0.80$ at every dose, for $\\delta>$ %g cm" % nb_, "region B is empty above this thickness (%sa)" % FIGN["fig6_falsification"]],
        ["Self-limiting enhancer (Gate G2)", "measured enhancer clearance not at least about 3-fold above the reference 0.5 h$^{-1}$, or recovery time above 2.5 h", "without it P(reversible) is %s with confinement (P5) against %s (P1)" % (pct(AR["P5_no_self_limiting_PE"]["p_safe"], 0), pct(p1["p_safe"], 0))],
        ["Matrix-only negative control", "matrix without enhancer raises $F$ more than about 3-fold over peptide alone", "the model gives %s-fold (M0: %s from the thin niche, %s from the shield); a larger gain in vivo would contradict the requirement for the enhancer [@Buckley2018]" % (f(AR["M0_matrix_only"]["median"] / AR["R0_peptide_alone"]["median"], 0), f(DC_[LUM]["R0 + thin niche (delta 0.30)"]["gain"], 1), f(DC_[LUM]["R0 + protection 0.70"]["gain"], 1))],
        ["Toxicity threshold (Gate G2)", "damage (viability or TEER loss) below the efficacious concentration, $C_{tox}/EC_{50}<$ 1", "P(success) of P1 falls to %s at a ratio of %s (Table 8)" % (pct(RV_["0.1"]["P1_dual_lock"]["p_success"], 0), f(RV_["0.1"]["C_tox_over_EC50_median"], 2))],
        ["Efflux module (Gate G3)", "efflux ratio below 2 for the peptide, or no reduction by the modulator", "the module has no effect for non-substrates (tested)"],
        ["Mixing and stagnant layer", "effective dispersion below about 100 times the molecular diffusivity", "the reference baseline of 1.2% cannot be reproduced in the spatial model below this (Table 10)"],
        ["Reference performance", "reference tablet at or below about 1% (semaglutide-like)", "P(success) of P1 falls to %s (Table 11)" % pct(A1["P1"]["p_success"], 0)],
        ["Absolute performance (Gate G4)", "median $F<4$% with P(reversible) $\\ge0.8$ in the in vivo study", "the 4 to 5% target is not met"]]
tbl(13, "Predefined failure conditions derived from the model. Each row states a measurable condition under which the corresponding hypothesis, or the architecture as a whole, would be considered falsified.",
    ["Hypothesis or gate", "Falsified if", "Reference in this paper"], fail, [2.0, 4.2, 2.8], 14)

# --------------------------------------------------------------------------------------------------------------------- 7 discussion
add(("h1", "7. Discussion"))
add(("p", "The analysis supports a narrow, conditional conclusion. If a matrix can confine the effective niche to about 0.3 cm and keep it mixed, and if the enhancer clears within tens of minutes, an architecture of this type would raise the modelled probability of reaching 4%% bioavailability with reversible exposure from %s to %s, "
          "at half the enhancer dose and without any assumed potency advantage. That is the only sense in which the proposal is supported, and the figure falls to %s if the reference tablet performs like oral semaglutide. Four aspects deserve emphasis."
          % (pct(r1["p_success"], 0), pct(p1["p_success"], 0), pct(A1["P1"]["p_success"], 0))))
add(("p", "*What is new, and what is not.* Permeation enhancers, P-glycoprotein modulators, mucoadhesive matrices and hydrophobic ion pairing are established. The contribution here is quantitative: a way to state when co-localization pays, how much of the benefit depends on each design variable, and where the architecture fails. "
          "The decomposition of the effect of design choices (large for the niche, small for the modulator) differs from the decomposition of the uncertainty of the outcome (shared among peptide biology, enhancer pharmacology and architecture), and both matter for experimental planning. "
          "The paired ablations also change the design: a slower release rate than the reference lowers success, so the matrix should release quickly and synchronously into the thin niche."))
add(("p", "*Biological plausibility.* The model reproduces the qualitative features that motivated it (a local and transient effect; dependence on volume and wash-out), but it fails one qualitative check (no optimum enhancer dose) and has an upper tail beyond any published value that the uptake constraints shorten but do not remove. "
          "The matrix-only control gains %s-fold over peptide alone in every formulation of the transport (original, stagnant layer, solubility limit, both), so it is not an artefact of the missing uptake limits; it is a consequence of the concentration effect of a thin niche (%s-fold) and of the protease shield (%s-fold). "
          "Because absorption of semaglutide requires SNAC [@Buckley2018], the matrix alone should not be expected to reach this level in vivo, and the first experiment therefore includes a matrix-only arm with a 3-fold acceptance threshold. "
          "The alternative transport structures did not remove the advantage, but only structures with substantial mixing reproduce the published baseline, so mixing is a second quantity the matrix must deliver, and the transport structure has not been independently validated (Section 6.2)."
          % (f(CN_[CON]["matrix_only_gain"], 1), f(DC_[CON]["R0 + thin niche (delta 0.30)"]["gain"], 1), f(DC_[CON]["R0 + protection 0.70"]["gain"], 1))))  # REV2
add(("p", "*Limitations.* All common priors and all arm settings are grade C or D; the dose, solubility and diffusion priors added for the uptake-constrained model are grade C or D as well. The model is linear in dose, lumps the niche, treats efflux as a divisor on permeability, and has no regional, motility or food effects. The reversibility criterion is a design target, not a safety limit, and the 4% target is borrowed from a published challenge, not derived from a clinical need. "
          "The calibration anchor is a single literature number; only two structural alternatives were tested; and no component of the work has been tested experimentally. No peptide-specific or independent dataset was found that tests the transport mechanism, and the reversibility indicator rests on assumed recovery kinetics, because damage-hours rarely bind in the base case. Chronic exposure and tolerability, the main safety questions for permeation enhancers, are not modelled."))
add(("p", "*Implications.* The most useful outputs are the failure boundaries of Table 13: they translate the model into thresholds that can be measured with a release test in a mucus-mimetic layer, a Caco-2 recovery assay and a bidirectional transport screen before any animal study. "
          "If the thresholds are missed, the architecture should be abandoned in its present form; if they are met, the priors can be replaced by measurements and the analysis repeated. The nominal design is not the best design in the model: a credibility-constrained optimization evaluated on independent draws gains %s and %s points of P(success) (original and uptake-constrained model) through faster, more synchronous release; "
          "this identifies a direction, not a validated performance."
          % (f(100 * OP_[LUM + "|credible"]["paired_opt_minus_nominal_test"][0], 1), f(100 * OP_[OPC + "|credible"]["paired_opt_minus_nominal_test"][0], 1))))

add(("h1", "8. Experimental validation roadmap"))
road = [["1", "Co-localization and mixing (G1)", "Biorelevant release in a small-volume apparatus with a mucus-mimetic layer: effective niche thickness, co-release of enhancer and peptide, stagnant-layer thickness; include a matrix-only arm.", "$\\delta\\le$ %g cm, release mismatch $m\\le0.2$, matrix-only gain below 3-fold, effective mixing consistent with a dispersion above 100 times molecular diffusivity" % nm_],
        ["2", "Reversible enhancement (G2)", "Caco-2 and mucus-secreting co-culture: $P_{app}$, TEER recovery, LDH release, enhancer clearance in brush-border or pancreatic preparations; library of about twelve amphiphiles.", "composite $E_{max}/EC_{50}$ in the upper half of the prior; TEER recovery within 2 h; clearance half-life of tens of minutes"],
        ["3", "Efflux (G3)", "Bidirectional transport in MDCK-MDR1 and Caco-2 with and without modulator; glycyl-sarcosine competition counter-screen for PepT1.", "efflux ratio $\\ge2$ and a reduction to $\\le1.3$ at a TPGS-compatible concentration"],
        ["4", "In vivo absolute bioavailability (G4)", "Large-animal study (beagle or minipig) against an intravenous reference with the reference formulation in parallel; short tolerability read-out.", "median $F\\ge4$% with reversible exposure; otherwise falsified (Table 13)"]]
tbl(14, "Prioritized validation studies and acceptance criteria derived from the model. Each study replaces grade-D priors by measurements, after which the analysis is repeated.", ["Step", "Question", "Experiment", "Acceptance criterion"], road, [0.5, 2.0, 4.0, 2.9], 14)

add(("h1", "9. Conclusions"))
add(("p", "A mechanistic, evidence-graded Monte Carlo analysis of a co-localized, reversible enhancer and efflux-modulator architecture for oral peptides gives the following conditional results. The architecture raises the modelled probability of reaching 4%% bioavailability with predicted reversible exposure from %s to %s at half the enhancer dose "
          "(%s against %s in the uptake-constrained model) and dominates the reference in the benefit-risk plane; the probability is %s to %s if the reference tablet performs like oral semaglutide. "
          "The benefit comes mainly from confinement and survives the stagnant-layer and solubility limits; it requires a self-limiting enhancer to remain reversible and a toxicity threshold at or above the efficacious concentration, and it is conditional on the matrix delivering a well-mixed niche of about 0.3 cm. "
          "Efflux modulation adds a smaller, substrate-dependent gain. The uncertainty of the prediction is shared among peptide biology, enhancer pharmacology and architecture; the nominal design is not optimal, and an out-of-sample test of a credibility-constrained optimum shows a further gain; the upper tail of the prediction is not credible. "
          "The transport structure has not been independently validated, and the study defines the measurements that would falsify it. It does not show that any formulation achieves 4 to 5%% bioavailability."
          % (pct(r1["p_success"], 0), pct(p1["p_success"], 0), pct(CN_[CON]["P1_dual_lock"]["p_success"], 0), pct(CN_[CON]["R1_reference_PE"]["p_success"], 0), pct(CN_[LUM]["abc_0.4_1.0"]["P1_success"], 0), pct(CN_[CON]["abc_0.4_1.0"]["P1_success"], 0))))

add(("h1", "CRediT authorship contribution statement"))
add(("p0", "**Leon Sandler:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Writing – review and editing, Visualization, Project administration."))
add(("h1", "Declaration of competing interest"))
add(("p0", "The author declares that he has no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper."))
add(("h1", "Funding"))
add(("p0", "This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors."))
add(("h1", "Data availability"))
add(("p0", "All code, tests, numerical results (JSON and CSV), parameter provenance tables and figure scripts are available at %s and archived on Zenodo (software, DOI %s; manuscript preprint, DOI %s). "
           "Fixed seeds regenerate every number; the extended analysis runs in about 15 minutes on one CPU core (Python 3.11.9, NumPy 2.4.6, SciPy 1.17.1, Matplotlib 3.11.1). The Supplementary Information accompanies this article."
           % (GITHUB, DOIS["code"], DOIS["paper"])))
add(("h1", "Declaration of generative AI and AI-assisted technologies in the writing process"))
add(("p0", "During the preparation of this work the author used Claude (Anthropic) to help develop and test the simulation code, to run the analyses and to draft and revise the text and figures. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the publication."))
add(("refs",))
