# -*- coding: utf-8 -*-
"""Revision: sections 5.8 (uptake constraints, matrix-only control), 5.9 (reversibility endpoints) and 5.10 (optimized design out of sample). Executed inside build_manuscript.py."""
CN_ = J("constrained")["variants"]; DC_ = J("decomp"); RV_ = J("revers")["ctox_ratio"]; OP_ = J("optim")
OPC = "uptake-constrained"
LUM, CON, LAYN, SOLN, DMG = "lumped (original)", "uptake-constrained (layer + solubility)", "stagnant layer", "solubility limit", "uptake-constrained, damage at efficacious concentration"

# ------------------------------------------------------------------------------------------------ 5.8
add(("h2", "5.8 Uptake constraints and the matrix-only control"))
add(("p", "The matrix-only control (Section 5.2) raised the median bioavailability %s-fold without any enhancer, which could be an artefact of a compartment model in which the uptake rate constant $P_{eff}/(\\delta\\,ER)$ grows without limit as the niche is thinned. "
          "Two physical limits were therefore added (Section 3.5): a stagnant layer in series with the epithelium and a solubility limit on the dissolved peptide, alone and together, each variant recalibrated to the 1.2%% anchor. "
          "Neither limit removes the gain: the matrix-only gain is %s with the stagnant layer, %s with the solubility limit and %s with both, against %s in the original model (%s). "
          "The solubility limit rarely binds for the assumed dose and solubility (the baseline permeability is unchanged, %s times, and the median of P1 falls only from %s to %s); the stagnant layer lowers the enhanced permeability but requires a baseline permeability %s times the lumped value, "
          "and the matrix-only control still reaches %s (peptide alone %s). The gain is therefore not caused by the missing uptake limits."
          % (f(CN_[LUM]["matrix_only_gain"], 1), f(CN_[LAYN]["matrix_only_gain"], 1), f(CN_[SOLN]["matrix_only_gain"], 1), f(CN_[CON]["matrix_only_gain"], 1), f(CN_[LUM]["matrix_only_gain"], 1), FIGN["fig11_mechanism_controls"],
             f(CN_[SOLN]["p0_over_base"], 2), pct(CN_[LUM]["P1_dual_lock"]["median"], 2), pct(CN_[SOLN]["P1_dual_lock"]["median"], 2), f(CN_[CON]["p0_over_base"], 1), pct(CN_[CON]["M0_matrix_only"]["median"], 2), pct(CN_[CON]["R0_peptide_alone"]["median"], 2))))
dl_, dc_ = DC_[LUM], DC_[CON]
add(("p", "The decomposition in %s shows where the gain comes from. Thinning the niche alone gives %s-fold, the protease shield alone %s-fold and the slower wash-out alone %s-fold (original model; uptake-constrained: %s, %s and %s). A matrix without the shield gains %s-fold and one without niche thinning %s-fold. "
          "The gain is thus mostly a concentration effect (the flux across the membrane is proportional to the concentration at the membrane, and the same dose dissolved in a smaller volume is more concentrated), plus the protease shield. "
          "The shield is not an external assumption: local buffering that protects against degradation is part of the reported mechanism of SNAC [@Buckley2018], so the matrix-only arm should be read as \"matrix with shield\", not as an inert carrier. "
          "What remains untested is whether a real matrix confines the dose into a smaller, mixed volume at all, which is the content of Gate G1."
          % (TABN_SEC.get("decomp", "Table 7"), f(dl_["R0 + thin niche (delta 0.30)"]["gain"], 1), f(dl_["R0 + protection 0.70"]["gain"], 1), f(dl_["R0 + slow wash-out (0.8)"]["gain"], 1),
             f(dc_["R0 + thin niche (delta 0.30)"]["gain"], 1), f(dc_["R0 + protection 0.70"]["gain"], 1), f(dc_["R0 + slow wash-out (0.8)"]["gain"], 1), f(dl_["M0 without protection"]["gain"], 1), f(dl_["M0 without niche thinning"]["gain"], 1))))
rows = []
for k in ("R0 peptide alone", "R0 + thin niche (delta 0.30)", "R0 + slow wash-out (0.8)", "R0 + protection 0.70", "R0 + niche + wash-out", "M0 without protection", "M0 without niche thinning", "M0 matrix only (all matrix parameters)"):
    rows.append([k.replace("delta 0.30", "$\\delta$ 0.30"), pct(dl_[k]["median"], 2), f(dl_[k]["gain"], 1), pct(dc_[k]["median"], 2), f(dc_[k]["gain"], 1)])
tbl(7, "Decomposition of the matrix-only gain over peptide alone (no enhancer, no modulator): median $F$ and gain, original lumped model and uptake-constrained model (both recalibrated; 20,000 shared draws).",
    ["Arm (no enhancer, no modulator)", "Lumped: median $F$", "Gain", "Constrained: median $F$", "Gain"], rows, [3.6, 1.3, 0.8, 1.5, 0.8], 14)
fig_fixed("fig11_mechanism_controls", "Mechanism controls under the original and constrained transport formulations. (a) Matrix-only gain (median bioavailability relative to peptide alone; dashed line, 3-fold acceptance threshold of Table 13). "
                                      "(b) Median bioavailability of the matrix-only arm and of peptide alone. (c) Volume-only effects: change in P(success) when the thin niche is added to R1 or removed from P1.", 6.4)
vo = lambda n: (100 * (CN_[n]["R1_plus_delta"]["p_success"] - CN_[n]["R1_reference_PE"]["p_success"]), 100 * (CN_[n]["P1_minus_delta"]["p_success"] - CN_[n]["P1_dual_lock"]["p_success"]))
add(("p", "The benefit of confinement itself is robust to the constraints. Adding the thin niche to R1 changes P(success) by %s points in the original model and %s points in the uptake-constrained model, and removing it from P1 by %s and %s points; "
          "in both formulations the thin niche is the largest single contributor and a self-limiting enhancer is still required, because adding the niche without it lowers the probability of reversible exposure (to %s, uptake-constrained)."
          % tuple(["%+.1f" % x for x in (vo(LUM)[0], vo(CON)[0], vo(LUM)[1], vo(CON)[1])] + [pct(CN_[CON]["R1_plus_delta"]["p_safe"], 0)])))

# ------------------------------------------------------------------------------------------------ 5.9
add(("h2", "5.9 Reversibility endpoints"))
r10, r03, r01 = RV_["1.0"], RV_["0.3"], RV_["0.1"]
gp = r10["P1_dual_lock"]["grid"]; rc = r10["P1_dual_lock"]["p_recov_le"]
add(("p", "In the base case the damage-hours rule seldom binds (median damage-hours %s h for P1 and %s h for R1 against the 0.5 h limit), so reversibility is decided by the assumed recovery time: the probability that $t_{clear}+t_{rev}\\le$ 2, 3, 4 and 6 h is %s, %s, %s and %s for P1. "
          "The predicted exposure reversibility is therefore a statement about the assumed recovery kinetics, and it should be tied to measurable endpoints: $t_{rev}$ to the recovery of transepithelial electrical resistance (TEER) after washout, and $C_{tox}$ to a viability endpoint such as lactate dehydrogenase release. "
          "Both are surrogates of epithelial integrity and not safety outcomes. Rat mucosae show a reduced TEER and minor epithelial damage at the concentrations at which SNAC and sodium caprate increase fluxes [@Twarog2021], which argues against a toxicity threshold far above the efficacious concentration."
          % (f(r10["P1_dual_lock"]["median_DH"], 4), f(r10["R1_reference_PE"]["median_DH"], 4), pct(rc["2.0"], 0), pct(rc["3.0"], 0), pct(rc["4.0"], 0), pct(rc["6.0"], 0))))
rows = []
for m_, lab in (("1.0", "3.2 (base)"), ("0.3", "0.96 (damage at efficacy)"), ("0.1", "0.32"), ("0.03", "0.096")):
    v = RV_[m_]; rows.append([lab] + ["%s / %s / %s" % (pct(v[a]["p_success"], 0), pct(v[a]["p_safe"], 0), f(v[a]["median_DH"], 3)) for a in ("R1_reference_PE", "P1_dual_lock", "P5_no_self_limiting_PE")])
tbl(8, "Sensitivity of the reversibility endpoint to the ratio of the damage concentration to the effective concentration, $C_{tox}/EC_{50}$ (medians). Cells: P(success) / P(predicted reversible) / median damage-hours (h).",
    ["$C_{tox}/EC_{50}$", "R1 reference", "P1 Dual-Lock", "P5 no self-limiting enhancer"], rows, [2.2, 2.3, 2.3, 2.5], 15)
add(("p", "When damage occurs at the efficacious concentration ($C_{tox}/EC_{50}=$ %s), P(success) is %s for P1 and %s for R1, so the advantage survives, and the reference and the non-self-limiting arm are hit hardest (P5 %s). "
          "When damage occurs below the efficacious concentration ($C_{tox}/EC_{50}=$ %s) the median damage-hours of P1 is %s h, %s of its draws exceed the 0.5 h limit and P(success) falls to %s (R1 %s); at a ratio of %s no arm succeeds. "
          "Tightening the damage-hours limit has a similar effect: with 4 h recovery, P1 P(success) is %s at 0.5 h, %s at 0.02 h and %s at 0.001 h. The Dual-Lock result is thus conditional on a toxicity threshold that is at or above the efficacious concentration, and Gate G2 must measure it."
          % (f(r03["C_tox_over_EC50_median"], 2), pct(r03["P1_dual_lock"]["p_success"], 0), pct(r03["R1_reference_PE"]["p_success"], 0), pct(r03["P5_no_self_limiting_PE"]["p_success"], 0),
             f(r01["C_tox_over_EC50_median"], 2), f(r01["P1_dual_lock"]["median_DH"], 2), pct(r01["P1_dual_lock"]["p_DH_gt_0p5"], 0), pct(r01["P1_dual_lock"]["p_success"], 0), pct(r01["R1_reference_PE"]["p_success"], 0), f(RV_["0.03"]["C_tox_over_EC50_median"], 2),
             pct(gp["0.5|4.0"], 0), pct(gp["0.02|4.0"], 0), pct(gp["0.001|4.0"], 0))))

# ------------------------------------------------------------------------------------------------ 5.10
add(("h2", "5.10 Optimized design evaluated out of sample"))
mk = [k for k in OP_ if "|" not in k]
add(("p", "The nominal design was not optimized: its release rate (1.2 h$^{-1}$) lowers success (Section 5.2). To separate the original design from an optimized one without presenting an in-sample optimum as performance, the release rate, enhancer dose and wash-out were optimized on 8,000 training draws "
          "(72 combinations) and the nominal design, the optimized design and the reference were evaluated on 20,000 independent test draws (%s, %s). "
          "The optimization is stable (optimism %s points in the original model and %s in the uptake-constrained model; the effect is too large for 72 candidates to overfit), but the unconstrained optimum lies on the boundary of the grid (release rate %g h$^{-1}$, the largest value) "
          "and produces a median bioavailability of %s with %s of draws above 10%%, which is outside the credible range of Section 6.2. It is therefore not a design recommendation."
          % (FIGN["fig12_optimized_vs_nominal"], TABN_SEC.get("optim", "Table 9"), f(100 * (OP_[LUM]["optimized_train"]["p_success"] - OP_[LUM]["optimized_test"]["p_success"]), 1), f(100 * (OP_[OPC]["optimized_train"]["p_success"] - OP_[OPC]["optimized_test"]["p_success"]), 1),
             OP_[LUM]["best"]["k_rel"], pct(OP_[LUM]["optimized_test"]["median"]), pct(OP_[LUM]["optimized_test"]["frac_gt_10"], 0))))
rows = []
for k in mk:
    v = OP_[k]; kk = k + "|credible"; c = OP_.get(kk)
    rows.append([k.replace(" (original)", ""), "nominal (0.5 mmol, 1.2 h$^{-1}$, 0.8 h$^{-1}$)", "%s (test)" % pct(v["nominal_test"]["p_success"]), pct(v["nominal_test"]["p_safe"], 0), pct(v["nominal_test"]["median"]), pct(v["nominal_test"]["frac_gt_10"], 0)])
    rows.append(["", "unconstrained optimum (%g mmol, %g h$^{-1}$, %g h$^{-1}$)" % (v["best"]["dose"], v["best"]["k_rel"], v["best"]["k_dil"]), "%s train / %s test" % (pct(v["optimized_train"]["p_success"]), pct(v["optimized_test"]["p_success"])), pct(v["optimized_test"]["p_safe"], 0), pct(v["optimized_test"]["median"]), pct(v["optimized_test"]["frac_gt_10"], 0)])
    if c:
        rows.append(["", "credible optimum (%g mmol, %g h$^{-1}$, %g h$^{-1}$)" % (c["best"]["dose"], c["best"]["k_rel"], c["best"]["k_dil"]), "%s train / %s test" % (pct(c["optimized_train"]["p_success"]), pct(c["optimized_test"]["p_success"])), pct(c["optimized_test"]["p_safe"], 0), pct(c["optimized_test"]["median"]), pct(c["optimized_test"]["frac_gt_10"], 0)])
    rows.append(["", "reference R1", "%s (test)" % pct(v["reference_test"]["p_success"]), pct(v["reference_test"]["p_safe"], 0), pct(v["reference_test"]["median"], 2), pct(v["reference_test"]["frac_gt_10"], 0)])
tbl(9, "Nominal and optimized designs: P(success) on the training draws (in-sample) and on 20,000 independent test draws, P(predicted reversible), median $F$ and fraction of draws above 10%. "
       "Design variables: enhancer dose, release rate $k_{rel}$ and wash-out $k_{dil}$. The credible optimum maximizes training P(success) subject to $k_{rel}\\le2.5$ h$^{-1}$ and at most 30% of draws above 10%.",
    ["Model", "Design (dose, release, wash-out)", "P(success)", "P(reversible)", "Median $F$", "Draws >10%"], rows, [1.6, 3.4, 1.9, 1.0, 0.9, 0.9], 13)
fig_fixed("fig12_optimized_vs_nominal", "Nominal design and optimized design in the original and uptake-constrained models: (a) P(success) with 95% Monte Carlo error, (b) probability of predicted reversible exposure, (c) median bioavailability (dashed line, 4%). "
                                        "Hatched bars: optimized design on the training draws (in-sample); solid bars: nominal and optimized design on independent test draws.", 6.4)
ck = [k + "|credible" for k in mk if k + "|credible" in OP_]
if ck:
    c0, c1 = OP_[ck[0]], OP_[ck[-1]]
    add(("p", "Restricting the search to a credible region gives a design with a release rate of %g h$^{-1}$ (original model) and %g h$^{-1}$ (uptake-constrained), whose out-of-sample P(success) is %s and %s, against %s and %s for the nominal design; the paired gain on the test draws is %s and %s points. "
              "These are model-conditional improvements that still rely on the priors; they identify the direction of a better design (faster, synchronous release) without validating its absolute performance."
              % (c0["best"]["k_rel"], c1["best"]["k_rel"], pct(c0["optimized_test"]["p_success"], 0), pct(c1["optimized_test"]["p_success"], 0), pct(c0["nominal_test"]["p_success"], 0), pct(c1["nominal_test"]["p_success"], 0),
                 f(100 * c0["paired_opt_minus_nominal_test"][0], 1), f(100 * c1["paired_opt_minus_nominal_test"][0], 1))))
