# -*- coding: utf-8 -*-
"""Part 3: front matter (title, authors, abstract, keywords) inserted at the start of the block list. Executed by build_manuscript.py."""
nb_ = FA["delta_max_for_B"]
ABSTRACT = [
    "Oral absolute bioavailability of 2 to 4 kDa peptides stays near 1%. We built a mechanistic Monte Carlo model of a co-localized, reversible permeation-enhancer and efflux-modulator architecture (Dual-Lock) with explicit mass balances and evidence-graded, correlated priors, and asked when it would raise bioavailability to at least 4% with predicted reversible epithelial exposure. "
    "Arms shared random draws, so that ablations and negative controls isolate each component. Global sensitivity, benefit-risk and falsification analyses, spatial and uptake-constrained transport models, reversibility-endpoint scenarios, observation-conditioned analysis and an out-of-sample evaluation of an optimized design tested robustness.",
    "With the reference calibrated to 1.2%%, Dual-Lock reached a median of %s at half the enhancer dose and met the target in %s of draws against %s for the reference. Confinement caused most of the benefit: a matrix-only control gained %s-fold over peptide alone, from thinning the niche (%s-fold) and the protease shield, and the gain survived stagnant-layer and solubility limits. "
    "In the uptake-constrained model success was %s (reference %s); it was %s to %s if the reference performs like oral semaglutide, and the upper tail stayed beyond published values. The advantage survived damage at the efficacious concentration but not below it. "
    "A credible release-rate optimum added %s to %s points out of sample. The transport structure is not independently validated; failure boundaries define the first experiments. Results are conditional on author-assigned priors."
    % (pct(p1["median"]), pct(p1["p_success"], 0), pct(r1["p_success"], 0), f(CN_[LUM]["matrix_only_gain"], 0), f(DC_[LUM]["R0 + thin niche (delta 0.30)"]["gain"], 1), pct(CN_[CON]["P1_dual_lock"]["p_success"], 0), pct(CN_[CON]["R1_reference_PE"]["p_success"], 0),
       pct(CN_[LUM]["abc_0.4_1.0"]["P1_success"], 0), pct(CN_[CON]["abc_0.4_1.0"]["P1_success"], 0), f(100 * OP_[LUM + "|credible"]["paired_opt_minus_nominal_test"][0], 0), f(100 * OP_[OPC + "|credible"]["paired_opt_minus_nominal_test"][0], 0)),
]
abs_words = len(" ".join(ABSTRACT).split())
print("abstract words:", abs_words)
assert abs_words <= 250, abs_words
HIGHLIGHTS = ["Oral peptide enhancers act in a micro-niche; a mechanistic model quantifies it",
              "Paired Monte Carlo arms isolate what co-localization adds to bioavailability",
              "Matrix-only gain survives uptake limits; it is thinning plus protease shield",
              "Optimized release is tested out of sample, not reported as in-sample optimum",
              "Failure boundaries give measurable pass-fail criteria for the first experiments"]
assert all(len(h) <= 85 for h in HIGHLIGHTS), [len(h) for h in HIGHLIGHTS]
KEYWORDS = "Oral peptide delivery; Permeation enhancer; Efflux modulation; Monte Carlo simulation; Global sensitivity analysis; Falsifiability"
FRONT = [("title", TITLE), ("centered", "Leon Sandler", 24), ("centered", "Independent Researcher, Northbrook, Illinois, USA", 20),
         ("centered", "Corresponding author: sandler.leon@gmail.com · ORCID 0009-0007-4584-808X", 20), ("blank",), ("h1", "Abstract")]
FRONT += [("p0", a) for a in ABSTRACT]
FRONT += [("p0", "**Keywords:** " + KEYWORDS), ("p0", "**Code and data:** %s (software, DOI %s); manuscript preprint, DOI %s." % (GITHUB, DOIS["code"], DOIS["paper"]))]
B[0:0] = FRONT
json.dump({"highlights": HIGHLIGHTS, "abstract": ABSTRACT, "keywords": KEYWORDS, "abstract_words": abs_words}, open(os.path.join(OUT, "frontmatter.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
