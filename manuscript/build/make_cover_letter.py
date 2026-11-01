# -*- coding: utf-8 -*-
"""Cover letter and Highlights for the International Journal of Pharmaceutics (numbers from results_ext/*.json, DOIs from dois.json).

    python make_cover_letter.py <results_dir> <out_dir>
"""
import datetime
import json
import os
import sys

from docx import Document
from docx.shared import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2])
J = lambda n: json.load(open(os.path.join(RES, "ext_%s.json" % n), encoding="utf-8"))
D = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))
FM = json.load(open(os.path.join(OUT, "frontmatter.json"), encoding="utf-8"))
GITHUB = "https://github.com/sandlerleon/dual-lock-oral-peptide-mc"
TITLE = "Mechanistic Monte Carlo assessment of co-localized permeation enhancement and efflux modulation for oral peptide delivery"
PA, AB, FA, ST = J("paired"), J("abc"), J("fals"), J("struct")
r1, p1 = PA["arms"]["R1_reference_PE"], PA["arms"]["P1_dual_lock"]
A1 = AB["R1 F in [0.4, 1.0] % (oral-semaglutide-like)"]
pc = lambda x, d=0: ("%." + str(d) + "f%%") % (100 * x)

d = Document()
st = d.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(11); st.paragraph_format.space_after = Pt(8)
for s in d.sections:
    s.left_margin = s.right_margin = Inches(1.0); s.top_margin = s.bottom_margin = Inches(0.9)


def para(t, bold=False, after=None):
    p = d.add_paragraph(); r = p.add_run(t); r.bold = bold
    if after is not None:
        p.paragraph_format.space_after = Pt(after)
    return p


para("Leon Sandler", True, 0)
for t in ("Independent Researcher", "Northbrook, Illinois, United States", "sandler.leon@gmail.com", "ORCID 0009-0007-4584-808X"):
    para(t, after=0)
para("")
para(datetime.date.today().strftime("%d %B %Y").lstrip("0"))
para("The Editors", after=0); para("International Journal of Pharmaceutics", after=0); para("Elsevier")
para("Dear Editors,")
para("I submit the manuscript “%s” for consideration as a research article in the International Journal of Pharmaceutics." % TITLE)
para("Oral peptide products reach an absolute bioavailability of only about 1%, with high variability. The mechanistic literature on the best-characterized enhancer tablet shows that the effect is local and transient: enhancer, peptide and protection act only where and when they coincide. "
     "The manuscript asks, quantitatively, when a dosage-form architecture that deliberately co-localizes a reversible permeation enhancer, a P-glycoprotein modulator and a protease shield in a thin niche would raise bioavailability above 4% with reversible epithelial exposure, and when it would fail. "
     "It is a computational design-space study, and I have written it as one: it contains no experimental data, every probability is conditional on author-assigned, evidence-graded priors, and the reference arm is calibrated to a published baseline, so it is not presented as validation.")
para("The study contributes a mechanistic model with every equation, unit and initial condition stated; paired Monte Carlo arms in which ablations and negative controls isolate co-localization, enhancer clearance, protection and efflux modulation; global variance-based sensitivity analysis with a decomposition of the output uncertainty; correlated priors; "
     "a benefit-risk analysis and falsification maps with predefined failure conditions; a structural-uncertainty analysis that replaces the lumped niche by stagnant-layer and spatial diffusion-reaction models; an observation-conditioned analysis; and comparison with published human bioavailability data that were not used for calibration.")
CN, DC, OP = J("constrained")["variants"], J("decomp"), J("optim")  # REV2
LU, CO = "lumped (original)", "uptake-constrained (layer + solubility)"
para("The principal findings are conditional and, in part, unfavorable to the concept. With the reference calibrated to 1.2%%, the Dual-Lock architecture reaches a median bioavailability of %s at half the enhancer dose and meets the target with predicted reversible exposure in %s of draws against %s for the reference. "
     "In response to the issues that a rigorous reviewer would raise, the revised paper adds: (i) physical uptake limits (a stagnant layer and a solubility limit) and a decomposition of the matrix-only control, which shows that its %s-fold gain is not removed by these limits and arises mainly from thinning the niche (%s-fold) and the protease shield; "
     "(ii) biologically constrained recalibration, in which P(success) is %s in the uptake-constrained model (reference %s) and falls to %s to %s if the reference tablet performs like oral semaglutide; the upper tail of the predicted bioavailability remains beyond any published value; "
     "(iii) reversibility endpoints, showing that the advantage persists if damage occurs at the efficacious concentration but not below it; and (iv) an optimized design evaluated on independent draws (a gain of %s to %s points for a credibility-constrained optimum). "
     "The transport structure has not been independently validated, and the paper turns the results into measurable failure boundaries (for example an effective niche thickness above %g cm) for the first experiments."
     % (pc(p1["median"], 1), pc(p1["p_success"]), pc(r1["p_success"]), "%.0f" % CN[LU]["matrix_only_gain"], "%.1f" % DC[LU]["R0 + thin niche (delta 0.30)"]["gain"], pc(CN[CO]["P1_dual_lock"]["p_success"]), pc(CN[CO]["R1_reference_PE"]["p_success"]),
        pc(CN[LU]["abc_0.4_1.0"]["P1_success"]), pc(CN[CO]["abc_0.4_1.0"]["P1_success"]), "%.0f" % (100 * OP[LU + "|credible"]["paired_opt_minus_nominal_test"][0] if LU + "|credible" in OP else 0),
        "%.0f" % (100 * OP["uptake-constrained|credible"]["paired_opt_minus_nominal_test"][0]), FA["delta_max_for_p30"]))
para("Code, tests, all numerical results and the Supplementary Information are openly available at %s and are archived at Zenodo (software, doi:%s; manuscript preprint, doi:%s). Fixed seeds regenerate every number. "
     "The original model was re-run from a clean directory before any change; its headline results reproduced exactly, and the one discrepancy (tie handling in a rank correlation, up to 0.03 in 18 sensitivity coefficients) is reported in the repository." % (GITHUB, D["code"], D["paper"]))
para("A generative AI assistant (Claude, Anthropic) was used to help develop and test the simulation code, to run the analyses and to draft and revise the manuscript text and figures, as declared in the manuscript. I specified the research questions, reviewed and edited all content and take full responsibility for it. "
     "The manuscript is original, has not been published elsewhere and is not under consideration at any other journal. I am the sole author, I have no competing interests, and the work received no specific funding. "
     "An earlier, non-peer-reviewed description of the concept was prepared as a proposal for an industry open-innovation challenge; the commercial content of that proposal is not part of this manuscript.")
para("Thank you for considering this submission.")
para("Sincerely,")
para("Leon Sandler", after=0)
para("Independent Researcher, Northbrook, Illinois, USA", after=0)
d.save(os.path.join(OUT, "Cover_Letter_IJP.docx")); print("wrote cover letter")

h = Document(); st = h.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(11)
h.add_paragraph().add_run("Highlights").bold = True
for t in FM["highlights"]:
    h.add_paragraph(t, style="List Bullet")
h.save(os.path.join(OUT, "DualLock_IJP_Highlights.docx")); print("wrote highlights")
