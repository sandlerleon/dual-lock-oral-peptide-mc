# -*- coding: utf-8 -*-
"""Supplementary Information for the Dual-Lock manuscript, built from results_ext/*.json and tables/*.csv.

    python build_si.py <results_dir> <figures_dir> <out_dir>
"""
import csv
import json
import os
import sys

import build_infra as BI
from refs_lib import LIB

RES = os.path.abspath(sys.argv[1]); FIGDIR = os.path.abspath(sys.argv[2]); OUT = os.path.abspath(sys.argv[3])
HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(RES, "ext_%s.json" % n), encoding="utf-8"))
DOIS = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))
GITHUB = "https://github.com/sandlerleon/dual-lock-oral-peptide-mc"
TAB = os.path.join(RES, "tables")
rd = lambda n: list(csv.reader(open(os.path.join(TAB, n + ".csv"), encoding="utf-8")))


def g(x, d=3):
    try:
        return ("%." + str(d) + "g") % float(x)
    except Exception:
        return x


B = []; add = B.append
add(("title", "Supplementary Information"))
add(("centered", "Mechanistic Monte Carlo assessment of co-localized permeation enhancement and efflux modulation for oral peptide delivery", 22))
add(("centered", "Leon Sandler, Independent Researcher, Northbrook, Illinois, USA", 20))
add(("centered", "Code and data: %s (software DOI %s; preprint DOI %s)" % (GITHUB, DOIS["code"], DOIS["paper"]), 18))
add(("blank",))
add(("p0", "This file contains the parameter provenance, the complete arm results, the full sensitivity indices, the structural variants, the convergence studies, the list of verification tests and the literature values used in the plausibility checks. "
           "All numbers are generated from the result files in the repository (results/ext_*.json, results/tables/*.csv)."))

add(("h1", "Table S1. Parameter priors and provenance"))
rows = [[r[0], r[2], r[3], r[4], r[5]] for r in rd("S1_priors_provenance")[1:]]
add(("table", "**Table S1.** Common priors with evidence grade and basis (grade A, direct measurement; B, close analogue; C, order-of-magnitude literature range; D, author assumption). No input reaches grade A. "
              "Eight of twelve common inputs are grade D and four are grade C.", ["Symbol", "Unit", "Grade", "Prior", "Basis"], rows, [0.7, 0.7, 0.5, 2.0, 4.4], 14))
rows = [[r[0]] + [g(x) for x in r[1:]] for r in rd("S1b_arms_definition")[1:]]
add(("table", "**Table S1b.** Arm definitions (medians of the formulation-controlled inputs; all grade D).", ["Arm", "PE (mmol)", "Modulator (mmol)", "$\\delta$ (cm)", "$k_{dil}$ (h$^{-1}$)", "$m$", "$Pr$", "$k_{dep}$ (h$^{-1}$)", "$k_{rel}$ (h$^{-1}$)"], rows,
     [3.0, 0.8, 0.9, 0.7, 0.8, 0.6, 0.6, 0.8, 0.8], 14))

add(("h1", "Table S2. Complete arm results (paired draws)"))
rows = [[r[0], "%.2f" % (100 * float(r[1])), "%.2f" % (100 * float(r[2])), "%.1f" % (100 * float(r[3])), "%.1f" % (100 * float(r[4])), "%.1f" % (100 * float(r[5])), "%.1f" % (100 * float(r[6])),
         "%+.1f" % (100 * float(r[11])), "%+.1f" % (100 * float(r[14]))] for r in rd("S2_arms_paired")[1:]]
add(("table", "**Table S2.** All arms of the paired design (20,000 shared draws). Percentages; $\\Delta$ is the paired change in P(success) relative to R1 or P1 in percentage points.",
     ["Arm", "Median $F$", "$F$ 5%", "$F$ 95%", "P($F\\ge4$%)", "P(reversible)", "P(success)", "$\\Delta$ vs R1", "$\\Delta$ vs P1"], rows, [3.2, 0.9, 0.8, 0.8, 0.9, 1.0, 0.9, 0.9, 0.9], 14))

add(("h1", "Table S3. Sobol sensitivity indices"))
S3 = rd("S3_sobol_indices")
for arm, lab in (("P1_dual_lock", "P1 Dual-Lock"), ("R1_reference_PE", "R1 reference")):
    rows = []
    for r in S3[1:]:
        if r[0] == arm and r[1] == "log10_F" and r[2] == "parameter":
            rows.append([r[3], "%.3f" % float(r[4]), "[%.3f, %.3f]" % (float(r[5]), float(r[6])), "%.3f" % float(r[7]), "[%.3f, %.3f]" % (float(r[8]), float(r[9]))])
    rows.sort(key=lambda x: -float(x[3]))
    cap = "**Table S3%s.** First-order ($S_1$) and total-order ($S_T$) indices for $\log_{10}F$, %s, with 95%% bootstrap intervals ($N=2^{13}$ base samples)." % ("a" if arm.startswith("P1") else "b", lab)
    add(("table", cap, ["Input", "$S_1$", "95% CI", "$S_T$", "95% CI"], rows, [1.4, 0.9, 1.7, 0.9, 1.7], 14))
rows = []
for r in S3[1:]:
    if r[2] == "group":
        rows.append([r[0][:2], r[1].replace("log10_", "log10 "), r[3], "%.3f" % float(r[4]), "%.3f" % float(r[7]), "[%.3f, %.3f]" % (float(r[8]), float(r[9]))])
add(("table", "**Table S3c.** Group indices (closed = group alone; total = group including interactions).", ["Arm", "Output", "Group", "Closed", "Total", "95% CI (total)"], rows, [0.6, 1.2, 2.6, 0.9, 0.9, 1.6], 14))

add(("h1", "Table S4. Structural variants"))
rows = []
for r in rd("S4_structural_variants")[1:]:
    rows.append([r[0]] + ([g(r[1]), g(r[2]), "%.2f" % (100 * float(r[3])), "%.2f" % (100 * float(r[4])), "%.1f" % (100 * float(r[5])), "%.1f" % (100 * float(r[6]))] if r[1] != "not calibratable" else ["not calibratable", "", "", "", "", ""]))
add(("table", "**Table S4.** All structural variants, including the variants in which the baseline permeability was *not* recalibrated (P0 kept at the lumped value). Median $F$ and P(success) in %; $P_0$ in cm h$^{-1}$.",
     ["Variant", "$P_0$", "$P_0$/base", "R1 median $F$", "P1 median $F$", "R1 P(success)", "P1 P(success)"], rows, [3.4, 0.9, 0.8, 1.0, 1.0, 1.0, 1.0], 13))

add(("h1", "Convergence and numerics"))
cv = J("converge")
rows = []
for arm in ("R1_reference_PE", "P1_dual_lock"):
    for d in cv["running_" + arm]:
        rows.append([arm[:2], "{:,}".format(d["n"]), "%.3f" % (100 * d["median_F"]), "%.2f" % (100 * d["p_success"]), "%.2f" % (100 * d["p_F_ge_4"])])
add(("table", "**Table S5.** Running estimates (one stream of $10^5$ draws).", ["Arm", "Draws", "Median $F$ (%)", "P(success) (%)", "P($F\\ge4$%) (%)"], rows, [0.8, 1.2, 1.4, 1.4, 1.4], 14))
fig_s1 = os.path.join(FIGDIR, "figS1_numerics.png")
add(("fig", fig_s1, "**Fig. S1.** Numerical convergence of the median bioavailability: time step (left, P1) and horizon (right).", 5.4))

add(("h1", "Table S6. Verification tests"))
tests = [("Original model (oral_peptide_mc.py --test)", "peptide mass balance closes to $10^{-9}$; $F\\in[0,1]$; determinism; calibration of R1 within 5%; enhancer raises $F$ over peptide alone; efflux module never lowers $F$ and has no effect on non-substrates; time-step convergence; rank correlation of identical vectors."),
         ("Sampler equivalence", "the uniform-variate sampler reproduces the original P1 success probability ($|z|<3$) and median $F$ (within 6%)."),
         ("Mass balance (lumped)", "closure error below $10^{-9}$ in every arm of the paired design (largest: $4\\times10^{-15}$)."),
         ("Sobol estimators", "the Saltelli and Jansen estimators recover the analytic first- and total-order indices of the Ishigami function (maximum error below 0.04) and give zero for inert inputs."),
         ("Gaussian copula", "realised Spearman correlation equals $(6/\\pi)\\arcsin(\\rho/2)$ within 0.01."),
         ("Spatial model, well-mixed limit", "for a very large dispersion multiplier the spatial model reproduces the lumped median $F$ within 3%."),
         ("Spatial model, mass balance", "closure error below $10^{-7}$ (observed $10^{-12}$)."),
         ("Spatial model, convergence", "doubling the mesh and halving the time step changes the median $F$ by less than 3%."),
         ("Spatial model, monotonicity", "slower mixing lowers $F$ for the same draws."),
         ("Stagnant layer", "never raises $F$ for the same draws; mass balance closes to $10^{-9}$."),
         ("Module logic with the new sampler", "no effect of the modulator if the peptide is never a substrate; the enhancer and modulator multiply the matrix-only median $F$ by a factor of at least 4."),
         ("Solubility limit and input modifiers", "the solubility limit never raises $F$ for the same draws and the mass balance closes to $10^{-9}$; lowering $C_{tox}$ never lowers damage-hours."),
         ("Determinism and helpers", "identical results for a fixed seed; the paired-difference helper returns the exact contrast.")]
add(("table", "**Table S6.** Automated verification tests (python oral_peptide_mc.py --test: 9 tests; python run_extended.py --test: 19 tests). They verify the software, not the biology.", ["Test", "What is checked"], [list(t) for t in tests], [2.2, 6.6], 15))

add(("h1", "Table S7. Literature values used in the plausibility checks"))
lit = [["Oral semaglutide (SNAC) absolute bioavailability", "0.8% under the recommended dosing conditions; within-subject variability 137%; falls with larger water volume and rises with longer post-dose fasting", "Overgaard et al. 2021, abstract", "[@Overgaard2021]"],
       ["Oral semaglutide, SNAC dose", "exposure highest with 300 mg SNAC among the doses tested (150 to 600 mg); exposure proportional to the semaglutide dose", "Granhall et al. 2019, abstract", "[@Granhall2019]"],
       ["SNAC mechanism", "transcellular; confined to the area next to the tablet; transient; requires co-formulation; local buffering protects against degradation", "Buckley et al. 2018, abstract", "[@Buckley2018]"],
       ["Enhancer-based oral peptides, typical bioavailability", "about 1%, low and variable", "Brayden et al. 2020, abstract", "[@Brayden2020]"],
       ["Engineered GLP-1 agonist with enhancers (dogs)", "bioavailability about 6% (enteric-coated tablets, sodium chenodeoxycholate and propyl gallate)", "Pechenov et al. 2021, abstract", "[@Pechenov2021]"],
       ["Adherent mucus gel thickness (rat)", "jejunum about 123 $\\mu$m; colon about 830 $\\mu$m; firmly adherent layer 20 to 154 $\\mu$m (small intestine 20 $\\mu$m or absent)", "Atuma et al. 2001, abstract", "[@Atuma2001]"]]
add(("table", "**Table S7.** Values quoted from the abstracts of the sources (full texts were not read). The values are used only to place the model output in context and are not calibration data, except for the 1 to 2% baseline of R1.",
     ["Quantity", "Value", "Source", "Ref."], lit, [2.2, 4.0, 1.7, 1.2], 14))

# REV2: tables of the revision
add(("h1", "Table S8. Uptake-constrained variants, decomposition, reversibility endpoints and optimization"))
rows = []
for r in rd("S8_constrained_variants")[1:]:
    if r[2] in ("R1_reference_PE", "P1_dual_lock", "M0_matrix_only", "R0_peptide_alone"):
        rows.append([r[0].replace("uptake-constrained, damage at efficacious concentration", "constrained + damage at efficacy"), g(r[1], 3), r[2], "%.2f" % (100 * float(r[3])), "%.1f" % (100 * float(r[4])), "%.1f" % (100 * float(r[6])), "%.1f" % (100 * float(r[8]))])
add(("table", "**Table S8a.** Model variants (each recalibrated to the 1.2% anchor): baseline permeability relative to the original, median and 95th percentile of $F$, P(success) and fraction of draws above 10%.",
     ["Variant", "$P_0$/base", "Arm", "Median $F$ (%)", "95th pct. (%)", "P(success) (%)", "$F>10$% (%)"], rows, [2.6, 0.8, 1.9, 0.9, 0.9, 0.9, 0.9], 13))
rows = [[r[0].replace("uptake-constrained, damage at efficacious concentration", "constrained + damage at efficacy"), r[1], "%.1f" % (100 * float(r[2])), "%.1f" % (100 * float(r[3])), "%.2f" % (100 * float(r[4]))] for r in rd("S8b_constrained_restrictions")[1:]]
add(("table", "**Table S8b.** Dual-Lock P1 under the restrictions of the biologically constrained recalibration.", ["Variant", "Restriction", "Draws retained (%)", "P1 P(success) (%)", "P1 median $F$ (%)"], rows, [2.6, 2.0, 1.1, 1.1, 1.1], 13))
rows = [[r[0], r[1], "%.2f" % (100 * float(r[2])), "%.2f" % (100 * float(r[3])), "%.1f" % float(r[4])] for r in rd("S9_matrix_only_decomposition")[1:]]
add(("table", "**Table S9.** Decomposition of the matrix-only gain (no enhancer, no modulator).", ["Model", "Arm", "Median $F$ (%)", "95th pct. (%)", "Gain"], rows, [2.2, 3.4, 1.0, 1.0, 0.7], 13))
rows = [[g(r[0], 3), r[1][:2], r[2], r[3], "%.1f" % (100 * float(r[4]))] for r in rd("S10_reversibility_grid")[1:] if r[1] == "P1_dual_lock" and r[3] in ("3.0", "4.0", "6.0")]
add(("table", "**Table S10.** P(success) of P1 for combinations of the damage-hours limit and the recovery limit, at several $C_{tox}/EC_{50}$ ratios (medians).", ["$C_{tox}/EC_{50}$", "Arm", "DH max (h)", "Recovery max (h)", "P(success) (%)"], rows, [1.4, 0.7, 1.2, 1.4, 1.2], 13))
rows = [[r[0], r[1], r[2][:60], "%.1f" % (100 * float(r[3])), "%.1f" % (100 * float(r[4])), "%.2f" % (100 * float(r[5])), "%.1f" % (100 * float(r[6])), "%.1f" % (100 * float(r[7]))] for r in rd("S11_optimization")[1:]]
add(("table", "**Table S11.** Nominal and optimized designs on training and test draws (design: release rate, dose, wash-out).", ["Model", "Set", "Design", "P(success)", "P(reversible)", "Median $F$", "95th pct.", "$F>10$%"], rows, [1.6, 1.4, 2.6, 0.9, 0.9, 0.8, 0.8, 0.8], 12))

add(("refs",))

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    cit = BI.Harvard(LIB)
    BI.scan_blocks(cit, B)
    # references used only in the SI table
    for k in list(cit.lib):
        pass
    unused = [k for k in LIB if k not in cit.used]
    for k in unused:
        cit.used.discard(k)
    keep = {}
    out = os.path.join(OUT, "DualLock_IJP_Supplementary_Information.docx")
    # restrict the reference list to the cited keys
    cit.lib = {k: v for k, v in LIB.items() if k in cit.used}
    BI.build(out, B, cit)
    print("wrote", out)
