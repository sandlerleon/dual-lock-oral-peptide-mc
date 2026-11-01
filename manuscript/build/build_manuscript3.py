# -*- coding: utf-8 -*-
"""Part 3: front matter (title, authors, abstract, keywords) inserted at the start of the block list. Executed by build_manuscript.py."""
nb_ = FA["delta_max_for_B"]
ABSTRACT = [
    "Oral absolute bioavailability of 2 to 4 kDa peptides stays near 1%. The best-characterized enhancer tablet acts only where and when enhancer, peptide and protection coincide, which suggests that the architecture of the dosage form is a design variable. "
    "We built a mechanistic Monte Carlo model of a co-localized, reversible permeation-enhancer and efflux-modulator architecture (Dual-Lock) with explicit mass balances, evidence-graded and correlated priors, and asked when it would raise bioavailability to at least 4% with reversible epithelial exposure. "
    "All arms shared the same random draws, so that ablations and negative controls isolate each component. Global sensitivity analysis, benefit-risk analysis, falsification maps, alternative transport models, observation-conditioned analysis and comparison with published human data tested robustness.",
    "With the reference calibrated to 1.2%%, Dual-Lock reached a median of %s at half the enhancer dose and met the target with reversible exposure in %s of draws against %s for the reference. Confinement accounted for most of the benefit (removing it cost %s of success) but required a self-limiting enhancer to remain reversible. "
    "Output uncertainty was shared among peptide biology, enhancer pharmacology and architecture. The advantage persisted in every transport structure that reproduced the baseline, but fell to %s if the reference performs like oral semaglutide, and the upper tail and the enhancer-dose response disagreed with published data. "
    "Failure boundaries (niche thickness above %g cm, an unmixed niche, a matrix-only gain above 3-fold) define the first experiments. The results are conditional on author-assigned priors and are not evidence that any formulation works."
    % (pct(p1["median"]), pct(p1["p_success"], 0), pct(r1["p_success"], 0), "%d points" % round(abs(100 * vp1["P1_minus_delta"]["d_success"][0])), pct(A1["P1"]["p_success"], 0), FA["delta_max_for_p30"]),
]
abs_words = len(" ".join(ABSTRACT).split())
print("abstract words:", abs_words)
assert abs_words <= 250, abs_words
HIGHLIGHTS = ["Oral peptide enhancers act in a micro-niche; a mechanistic model quantifies it",
              "Paired Monte Carlo arms isolate what co-localization adds to bioavailability",
              "Variance-based sensitivity ranks niche thickness and enhancer potency first",
              "Spatial and stagnant-layer models test whether the predicted advantage survives",
              "Failure boundaries give measurable pass-fail criteria for the first experiments"]
assert all(len(h) <= 85 for h in HIGHLIGHTS), [len(h) for h in HIGHLIGHTS]
KEYWORDS = "Oral peptide delivery; Permeation enhancer; Efflux modulation; Monte Carlo simulation; Global sensitivity analysis; Falsifiability"
FRONT = [("title", TITLE), ("centered", "Leon Sandler", 24), ("centered", "Independent Researcher, Northbrook, Illinois, USA", 20),
         ("centered", "Corresponding author: sandler.leon@gmail.com · ORCID 0009-0007-4584-808X", 20), ("blank",), ("h1", "Abstract")]
FRONT += [("p0", a) for a in ABSTRACT]
FRONT += [("p0", "**Keywords:** " + KEYWORDS), ("p0", "**Code and data:** %s (software, DOI %s); manuscript preprint, DOI %s." % (GITHUB, DOIS["code"], DOIS["paper"]))]
B[0:0] = FRONT
json.dump({"highlights": HIGHLIGHTS, "abstract": ABSTRACT, "keywords": KEYWORDS, "abstract_words": abs_words}, open(os.path.join(OUT, "frontmatter.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
