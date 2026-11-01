# -*- coding: utf-8 -*-
"""Revision: biologically constrained recalibration and the status of independent validation (inserted in section 6.2). Executed inside build_manuscript.py."""
rows = []
for k in (LUM, CON, DMG):
    v = CN_[k]
    rows.append([k.replace(" (original)", "").replace("uptake-constrained, damage at efficacious concentration", "constrained + damage at efficacy").replace("uptake-constrained (layer + solubility)", "uptake-constrained"),
                 "%sx" % f(v["p0_over_base"], 1), pct(v["P1_dual_lock"]["median"]), pct(v["P1_dual_lock"]["p95"], 0), pct(v["P1_dual_lock"]["frac_gt_6"], 0), pct(v["P1_dual_lock"]["p_success"], 0), pct(v["abc_0.4_2.0"]["P1_success"], 0),
                 pct(v["abc_0.4_1.0"]["P1_success"], 0), pct(v["tail_rejected_F_le_10"]["P1_success"], 0), pct(v["abc_and_tail"]["P1_success"], 0)])
tbl(12, "Biologically constrained recalibration of Dual-Lock P1: model variant, recalibrated baseline permeability relative to the original, median and 95th percentile of $F$, fraction of draws above 6%, and P(success) for all draws and after each restriction "
        "(R1 in 0.4 to 2.0%; R1 in 0.4 to 1.0%; no arm above 10%; R1 in 0.4 to 2.0% and no arm above 10%).",
    ["Model", "$P_0$", "Median $F$", "95th pct.", "$F>6$%", "All draws", "R1 0.4-2.0%", "R1 0.4-1.0%", "$F\\le10$%", "both"], rows, [2.4, 0.7, 0.9, 0.8, 0.8, 0.9, 1.0, 1.0, 0.9, 0.8], 13)
fig_fixed("fig10_constrained", "Biologically constrained recalibration. (a) P(success) of R1 and P1 in the original model, with each physical limit, with both limits, and with damage at the efficacious concentration. "
                               "(b) P(success) of P1 under the restrictions of Table 12 in the original and in the uptake-constrained model.", 6.4)
L1, C1, C2 = CN_[LUM], CN_[CON], CN_[DMG]
add(("p", "A peptide-specific calibration is not possible, because no measurements exist for the peptide that the architecture would carry; the nearest peptide-specific anchor is oral semaglutide (0.8%%, with a within-subject variability of 137%%) [@Overgaard2021]. "
          "Table 12 therefore reports how much of the advantage survives progressively stronger constraints. In the uptake-constrained model P1 has a median of %s and a 95th percentile of %s (original model %s and %s), and the fraction of draws above 6%% is %s (original %s): the physical limits shorten the tail but do not remove it. "
          "If values above 10%% are treated as non-credible and removed (%s of the draws are retained), P(success) is %s in the original and %s in the uptake-constrained model; if, in addition, the reference tablet must lie in 0.4 to 2.0%%, it is %s and %s; "
          "if the reference is oral semaglutide-like (0.4 to 1.0%%), %s and %s; and with damage at the efficacious concentration, %s. Between roughly one fifth and one half of the draws therefore still meet the target under every constraint, but the advantage over the reference is "
          "smaller and the absolute figures are lower than the headline value."
          % (pct(C1["P1_dual_lock"]["median"]), pct(C1["P1_dual_lock"]["p95"], 0), pct(L1["P1_dual_lock"]["median"]), pct(L1["P1_dual_lock"]["p95"], 0), pct(C1["P1_dual_lock"]["frac_gt_6"], 0), pct(L1["P1_dual_lock"]["frac_gt_6"], 0),
             pct(C1["tail_rejected_F_le_10"]["frac"], 0), pct(L1["tail_rejected_F_le_10"]["P1_success"], 0), pct(C1["tail_rejected_F_le_10"]["P1_success"], 0), pct(L1["abc_and_tail"]["P1_success"], 0), pct(C1["abc_and_tail"]["P1_success"], 0),
             pct(L1["abc_0.4_1.0"]["P1_success"], 0), pct(C1["abc_0.4_1.0"]["P1_success"], 0), pct(C2["P1_dual_lock"]["p_success"], 0))))
add(("p", "*Status of independent validation.* The transport structure has not been validated independently. The published measurements found in the sources read do not contain the variables needed (peptide permeability as a function of enhancer concentration and exposure time through a mucus-like layer) in a form that can be compared quantitatively. "
          "Three consistency checks that did not set any parameter are available. (i) In isolated rat jejunal and colonic mucosae, 20 mM SNAC was needed to reach the efficacy of 10 mM sodium caprate [@Twarog2021], consistent in order of magnitude with the millimolar-to-tens-of-millimolar $EC_{50}$ prior (median 25 mM). "
          "(ii) Measured adherent mucus thicknesses of 20 to 830 $\\mu$m [@Atuma2001] bracket the stagnant-layer prior. (iii) The weak binding of exenatide to SNAC and sodium caprate (dissociation constants of about 10 to 100 $\\mu$M) [@Twarog2022] supports omitting complexation. "
          "None of these tests the Dual-Lock transport mechanism. The prespecified test is the mucus-mimetic release experiment of Gate G1 (Section 8), and until it is done the transport structure remains an assumption, which is why the conclusions are stated as conditional."))
