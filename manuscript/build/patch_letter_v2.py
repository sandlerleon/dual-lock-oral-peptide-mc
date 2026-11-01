# -*- coding: utf-8 -*-
import os
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_cover_letter.py")
s = open(p, encoding="utf-8").read()
if "REV2" in s:
    raise SystemExit("already")
i = s.index('para("The principal findings'); j = s.index('para("Code, tests')
new = '''CN, DC, OP = J("constrained")["variants"], J("decomp"), J("optim")  # REV2
LU, CO = "lumped (original)", "uptake-constrained (layer + solubility)"
para("The principal findings are conditional and, in part, unfavorable to the concept. With the reference calibrated to 1.2%%, the Dual-Lock architecture reaches a median bioavailability of %s at half the enhancer dose and meets the target with predicted reversible exposure in %s of draws against %s for the reference. "
     "In response to the issues that a rigorous reviewer would raise, the revised paper adds: (i) physical uptake limits (a stagnant layer and a solubility limit) and a decomposition of the matrix-only control, which shows that its %s-fold gain is not removed by these limits and arises mainly from thinning the niche (%s-fold) and the protease shield; "
     "(ii) biologically constrained recalibration, in which P(success) is %s in the uptake-constrained model (reference %s) and falls to %s to %s if the reference tablet performs like oral semaglutide; the upper tail of the predicted bioavailability remains beyond any published value; "
     "(iii) reversibility endpoints, showing that the advantage persists if damage occurs at the efficacious concentration but not below it; and (iv) an optimized design evaluated on independent draws (a gain of %s to %s points for a credibility-constrained optimum). "
     "The transport structure has not been independently validated, and the paper turns the results into measurable failure boundaries (for example an effective niche thickness above %g cm) for the first experiments."
     % (pc(p1["median"], 1), pc(p1["p_success"]), pc(r1["p_success"]), "%.0f" % CN[LU]["matrix_only_gain"], "%.1f" % DC[LU]["R0 + thin niche (delta 0.30)"]["gain"], pc(CN[CO]["P1_dual_lock"]["p_success"]), pc(CN[CO]["R1_reference_PE"]["p_success"]),
        pc(CN[LU]["abc_0.4_1.0"]["P1_success"]), pc(CN[CO]["abc_0.4_1.0"]["P1_success"]), "%.0f" % (100 * OP[LU + "|credible"]["paired_opt_minus_nominal_test"][0] if LU + "|credible" in OP else 0),
        "%.0f" % (100 * OP["uptake-constrained|credible"]["paired_opt_minus_nominal_test"][0]), FA["delta_max_for_p30"]))
'''
s = s[:i] + new + s[j:]
s = s.replace("meets the target with reversible exposure", "meets the target with predicted reversible exposure")
open(p, "w", encoding="utf-8").write(s)
print("ok")
