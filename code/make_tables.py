# -*- coding: utf-8 -*-
"""make_tables.py -- CSV tables for the Supplementary Information from results_ext/*.json (parameter provenance, arms, Sobol, structural variants, grids)."""
import csv
import json
import os
import sys

import oral_peptide_mc as M

OUT = sys.argv[1] if len(sys.argv) > 1 else "results_ext"
TAB = os.path.join(OUT, "tables"); os.makedirs(TAB, exist_ok=True)
J = lambda n: json.load(open(os.path.join(OUT, "ext_%s.json" % n)))


def w(name, header, rows):
    with open(os.path.join(TAB, name + ".csv"), "w", newline="", encoding="utf-8") as fh:
        cw = csv.writer(fh); cw.writerow(header); cw.writerows(rows)
    print("wrote", name, len(rows))


BASIS = {"P0": "Calibrated so that the median of the reference arm R1 equals 1.2 % (mid-point of the 1-2 % enhancer-tablet range, Brayden et al. 2020; oral semaglutide 0.8 %, Overgaard et al. 2021). Dispersion: author assumption.",
         "k_deg": "Author assumption (no direct measurement for an enhancer niche); order of 1 h^-1 for luminal proteolysis of an unprotected peptide.",
         "S_post": "Author assumption (intra-enterocyte handling and first pass).",
         "ER0": "Caco-2 bidirectional efflux ratios for P-glycoprotein substrates, order of magnitude (Collnot et al. 2006); applies only if the peptide is a substrate.",
         "p_sub": "Scenario axis, not a prior (peptide-specific; unknown for any given 2-4 kDa peptide).",
         "A_m": "Author assumption: mucosal contact area of the tablet niche.",
         "EC50": "Millimolar range of SNAC and sodium caprate effective concentrations (Maher et al. 2016; Twarog et al. 2019).",
         "Emax": "Author assumption (maximal fold-increase of peptide permeability).",
         "h": "Author assumption (Hill coefficient).",
         "C_tox": "Author assumption (half-damage concentration of the enhancer).",
         "t_rev": "Author assumption (epithelial recovery time after enhancer clearance).",
         "IC50": "Piperine 15-74 uM in Caco-2 (Bhardwaj et al. 2002); TPGS sub-millimolar (Collnot et al. 2006); wide dispersion."}
UNIT = {k: v["unit"] for k, v in M.PRIORS.items()}
rows = []
for k, v in M.PRIORS.items():
    if v["kind"] == "ln":
        pr = "lognormal median %s, sigma_ln %s" % ("0.009746 (calibrated)" if k == "P0" else v["med"], v["sig"])
    elif v["kind"] == "beta":
        pr = "beta mean %s, sd %s" % (v["mean"], v["sd"])
    else:
        pr = "scenario axis, base %s" % v["val"]
    rows.append([k, v["note"], UNIT[k], v["grade"], pr, BASIS[k]])
w("S1_priors_provenance", ["symbol", "meaning", "unit", "evidence_grade", "prior", "basis"], rows)
w("S1b_arms_definition", ["arm", "pe_mmol", "inh_mmol", "delta_cm", "k_dil_per_h", "mism_sigma_ln", "protection", "k_dep_per_h", "k_rel_per_h"],
  [[k, a["pe"], a["inh"], a["delta"], a["k_dil"], a["mism"], a["prot"], a["k_dep"], a["k_rel"]] for k, a in M.ARMS.items()])
P = J("paired")
w("S2_arms_paired", ["arm", "median_F", "p05", "p95", "P_F_ge_4", "P_safe", "P_success", "se_success", "median_DH", "p90_DH", "closure_error", "d_success_vs_R1", "d_success_vs_R1_lo", "d_success_vs_R1_hi", "d_success_vs_P1", "d_success_vs_P1_lo", "d_success_vs_P1_hi"],
  [[k, s["median"], s["p05"], s["p95"], s["p_F_ge_4"], s["p_safe"], s["p_success"], s["se"], s["median_DH"], s["p90_DH"], s["closure_err"], *P["paired_vs_R1"][k]["d_success"], *P["paired_vs_P1"][k]["d_success"]] for k, s in P["arms"].items()])
S = J("sobol"); rows = []
for arm in ("P1_dual_lock", "R1_reference_PE"):
    for o in ("log10_F", "success", "log10_DH"):
        for k, v in S[arm][o]["params"].items():
            rows.append([arm, o, "parameter", k, v["S1"], v["S1_ci"][0], v["S1_ci"][1], v["ST"], v["ST_ci"][0], v["ST_ci"][1]])
        for k, v in S[arm][o]["groups"].items():
            rows.append([arm, o, "group", k, v["closed"], v["closed_ci"][0], v["closed_ci"][1], v["total"], v["total_ci"][0], v["total_ci"][1]])
w("S3_sobol_indices", ["arm", "output", "level", "name", "S1_or_closed", "lo", "hi", "ST_or_total", "lo_T", "hi_T"], rows)
St = J("struct")["variants"]; rows = []
for k, v in St.items():
    if v.get("calibratable", True) is False:
        rows.append([k, "not calibratable", "", "", "", "", ""]); continue
    rows.append([k, v["p0"], v["p0_over_base"], v["R1_reference_PE"]["median"], v["P1_dual_lock"]["median"], v["R1_reference_PE"]["p_success"], v["P1_dual_lock"]["p_success"]])
w("S4_structural_variants", ["variant", "P0_cm_per_h", "P0_over_base", "R1_median_F", "P1_median_F", "R1_P_success", "P1_P_success"], rows)
Fa = J("fals")
for name, key, rl, cl in (("S5a_grid_delta_dose_P_success", "grid_delta_dose", Fa["doses"], Fa["deltas"]), ("S5b_grid_delta_kdil_P_success", "grid_delta_kdil", Fa["k_dils"], Fa["deltas"])):
    G = Fa[key]["p_success"]
    w(name, ["row\\delta_cm"] + cl, [[r] + G[i] for i, r in enumerate(rl)])
Pa = J("pareto"); rows = []
for fam, lst in Pa["families"].items():
    for s in lst:
        rows.append([fam, s["dose"], s["median"], s["p_F_ge_4"], s["p_safe"], s["p_success"], s["median_DH"], s["p90_DH"]])
w("S6_pareto_dose_sweeps", ["family", "dose_mmol", "median_F", "P_F_ge_4", "P_safe", "P_success", "median_DH", "p90_DH"], rows)
