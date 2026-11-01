# -*- coding: utf-8 -*-
"""Revision text patches for build_manuscript2.py / build_manuscript3.py (discussion, limitations, conclusions, abstract)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
p2 = os.path.join(HERE, "build_manuscript2.py"); p3 = os.path.join(HERE, "build_manuscript3.py")
s = open(p2, encoding="utf-8").read()
if "REV2" in s:
    raise SystemExit("already patched")


def rep(s, a, b):
    assert a in s, a[:70]
    return s.replace(a, b, 1)


# --- discussion: biological plausibility paragraph
i = s.index('add(("p", "*Biological plausibility.*'); j = s.index('add(("p", "*Limitations.*')
s = s[:i] + '''add(("p", "*Biological plausibility.* The model reproduces the qualitative features that motivated it (a local and transient effect; dependence on volume and wash-out), but it fails one qualitative check (no optimum enhancer dose) and has an upper tail beyond any published value that the uptake constraints shorten but do not remove. "
          "The matrix-only control gains %s-fold over peptide alone in every formulation of the transport (original, stagnant layer, solubility limit, both), so it is not an artefact of the missing uptake limits; it is a consequence of the concentration effect of a thin niche (%s-fold) and of the protease shield (%s-fold). "
          "Because absorption of semaglutide requires SNAC [@Buckley2018], the matrix alone should not be expected to reach this level in vivo, and the first experiment therefore includes a matrix-only arm with a 3-fold acceptance threshold. "
          "The alternative transport structures did not remove the advantage, but only structures with substantial mixing reproduce the published baseline, so mixing is a second quantity the matrix must deliver, and the transport structure has not been independently validated (Section 6.2)."
          % (f(CN_[CON]["matrix_only_gain"], 1), f(DC_[CON]["R0 + thin niche (delta 0.30)"]["gain"], 1), f(DC_[CON]["R0 + protection 0.70"]["gain"], 1))))  # REV2
''' + s[j:]

# --- limitations
s = rep(s, 'All twelve common priors and all arm settings are grade C or D.', 'All common priors and all arm settings are grade C or D; the dose, solubility and diffusion priors added for the uptake-constrained model are grade C or D as well.')
s = rep(s, "and no component of the work has been tested experimentally.", "and no component of the work has been tested experimentally. No peptide-specific or independent dataset was found that tests the transport mechanism, and the reversibility indicator rests on assumed recovery kinetics, because damage-hours rarely bind in the base case.")

# --- implications: optimization
i = s.index('add(("p", "*Implications.*'); j = s.index('add(("h1", "8. Experimental validation roadmap"))')
s = s[:i] + '''add(("p", "*Implications.* The most useful outputs are the failure boundaries of Table 13: they translate the model into thresholds that can be measured with a release test in a mucus-mimetic layer, a Caco-2 recovery assay and a bidirectional transport screen before any animal study. "
          "If the thresholds are missed, the architecture should be abandoned in its present form; if they are met, the priors can be replaced by measurements and the analysis repeated. The nominal design is not the best design in the model: a credibility-constrained optimization evaluated on independent draws gains %s and %s points of P(success) (original and uptake-constrained model) through faster, more synchronous release; "
          "this identifies a direction, not a validated performance."
          % (f(100 * OP_[LUM + "|credible"]["paired_opt_minus_nominal_test"][0], 1), f(100 * OP_[CON + "|credible"]["paired_opt_minus_nominal_test"][0], 1))))

''' + s[j:]

# --- failure table: matrix-only row and reversibility row
s = rep(s, '"the model gives %s-fold (M0); a larger gain in vivo would contradict the requirement for the enhancer [@Buckley2018]" % f(AR["M0_matrix_only"]["median"] / AR["R0_peptide_alone"]["median"], 0)]',
        '"the model gives %s-fold (M0: %s from the thin niche, %s from the shield); a larger gain in vivo would contradict the requirement for the enhancer [@Buckley2018]" % (f(AR["M0_matrix_only"]["median"] / AR["R0_peptide_alone"]["median"], 0), f(DC_[LUM]["R0 + thin niche (delta 0.30)"]["gain"], 1), f(DC_[LUM]["R0 + protection 0.70"]["gain"], 1))]')
s = rep(s, '        ["Efflux module (Gate G3)"', '        ["Toxicity threshold (Gate G2)", "damage (viability or TEER loss) below the efficacious concentration, $C_{tox}/EC_{50}<$ 1", "P(success) of P1 falls to %s at a ratio of %s (Table 8)" % (pct(RV_["0.1"]["P1_dual_lock"]["p_success"], 0), f(RV_["0.1"]["C_tox_over_EC50_median"], 2))],\n        ["Efflux module (Gate G3)"')

# --- conclusions paragraph
i = s.index('add(("h1", "9. Conclusions"))'); j = s.index('add(("h1", "CRediT authorship contribution statement"))')
s = s[:i] + '''add(("h1", "9. Conclusions"))
add(("p", "A mechanistic, evidence-graded Monte Carlo analysis of a co-localized, reversible enhancer and efflux-modulator architecture for oral peptides gives the following conditional results. The architecture raises the modelled probability of reaching 4%% bioavailability with predicted reversible exposure from %s to %s at half the enhancer dose "
          "(%s against %s in the uptake-constrained model) and dominates the reference in the benefit-risk plane; the probability is %s to %s if the reference tablet performs like oral semaglutide. "
          "The benefit comes mainly from confinement and survives the stagnant-layer and solubility limits; it requires a self-limiting enhancer to remain reversible and a toxicity threshold at or above the efficacious concentration, and it is conditional on the matrix delivering a well-mixed niche of about 0.3 cm. "
          "Efflux modulation adds a smaller, substrate-dependent gain. The uncertainty of the prediction is shared among peptide biology, enhancer pharmacology and architecture; the nominal design is not optimal, and an out-of-sample test of a credibility-constrained optimum shows a further gain; the upper tail of the prediction is not credible. "
          "The transport structure has not been independently validated, and the study defines the measurements that would falsify it. It does not show that any formulation achieves 4 to 5%% bioavailability."
          % (pct(r1["p_success"], 0), pct(p1["p_success"], 0), pct(CN_[CON]["P1_dual_lock"]["p_success"], 0), pct(CN_[CON]["R1_reference_PE"]["p_success"], 0), pct(CN_[LUM]["abc_0.4_1.0"]["P1_success"], 0), pct(CN_[CON]["abc_0.4_1.0"]["P1_success"], 0))))

''' + s[j:]
open(p2, "w", encoding="utf-8").write(s)

# --- abstract (part 3)
t = open(p3, encoding="utf-8").read()
i = t.index("ABSTRACT = ["); j = t.index("abs_words = ")
t = t[:i] + '''ABSTRACT = [
    "Oral absolute bioavailability of 2 to 4 kDa peptides stays near 1%. We built a mechanistic Monte Carlo model of a co-localized, reversible permeation-enhancer and efflux-modulator architecture (Dual-Lock) with explicit mass balances and evidence-graded, correlated priors, and asked when it would raise bioavailability to at least 4% with predicted reversible epithelial exposure. "
    "Arms shared random draws, so that ablations and negative controls isolate each component. Global sensitivity, benefit-risk and falsification analyses, spatial and uptake-constrained transport models, reversibility-endpoint scenarios, observation-conditioned analysis and an out-of-sample evaluation of an optimized design tested robustness.",
    "With the reference calibrated to 1.2%%, Dual-Lock reached a median of %s at half the enhancer dose and met the target in %s of draws against %s for the reference. Confinement caused most of the benefit: a matrix-only control gained %s-fold over peptide alone, from thinning the niche (%s-fold) and the protease shield, and the gain survived stagnant-layer and solubility limits. "
    "In the uptake-constrained model success was %s (reference %s); it was %s to %s if the reference performs like oral semaglutide, and the upper tail stayed beyond published values. The advantage survived damage at the efficacious concentration but not below it. "
    "A credible release-rate optimum added %s to %s points out of sample. The transport structure is not independently validated; failure boundaries define the first experiments. Results are conditional on author-assigned priors."
    % (pct(p1["median"]), pct(p1["p_success"], 0), pct(r1["p_success"], 0), f(CN_[LUM]["matrix_only_gain"], 0), f(DC_[LUM]["R0 + thin niche (delta 0.30)"]["gain"], 1), pct(CN_[CON]["P1_dual_lock"]["p_success"], 0), pct(CN_[CON]["R1_reference_PE"]["p_success"], 0),
       pct(CN_[LUM]["abc_0.4_1.0"]["P1_success"], 0), pct(CN_[CON]["abc_0.4_1.0"]["P1_success"], 0), f(100 * OP_[LUM + "|credible"]["paired_opt_minus_nominal_test"][0], 0), f(100 * OP_[CON + "|credible"]["paired_opt_minus_nominal_test"][0], 0)),
]
''' + t[j:]
t = t.replace('"Variance-based sensitivity ranks niche thickness and enhancer potency first"', '"Matrix-only gain survives uptake limits; it is thinning plus protease shield"')
t = t.replace('"Spatial and stagnant-layer models test whether the predicted advantage survives"', '"Optimized release is tested out of sample, not reported as in-sample optimum"')
open(p3, "w", encoding="utf-8").write(t)
print("patched")
