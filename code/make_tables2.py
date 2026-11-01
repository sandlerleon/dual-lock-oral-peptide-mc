# -*- coding: utf-8 -*-
"""Additional CSV tables for the revision (constrained variants, decomposition, reversibility grid, optimization). Run after make_tables.py."""
import csv
import json
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "results_ext"
TAB = os.path.join(OUT, "tables"); os.makedirs(TAB, exist_ok=True)
J = lambda n: json.load(open(os.path.join(OUT, "ext_%s.json" % n)))


def w(name, header, rows):
    with open(os.path.join(TAB, name + ".csv"), "w", newline="", encoding="utf-8") as fh:
        cw = csv.writer(fh); cw.writerow(header); cw.writerows(rows)
    print("wrote", name, len(rows))


C = J("constrained")["variants"]; rows = []
for v, d in C.items():
    for arm in ("R0_peptide_alone", "M0_matrix_only", "R1_reference_PE", "P1_dual_lock", "P1_minus_delta", "R1_plus_delta", "P3_no_colocalisation", "P5_no_self_limiting_PE", "R1_plus_confinement_bundle"):
        s = d[arm]; rows.append([v, "%.3f" % d["p0_over_base"], arm, s["median"], s["p95"], s["p_safe"], s["p_success"], s["frac_gt_6"], s["frac_gt_10"]])
w("S8_constrained_variants", ["variant", "P0_over_base", "arm", "median_F", "p95_F", "P_safe", "P_success", "frac_F_gt_6", "frac_F_gt_10"], rows)
rows = []
for v, d in C.items():
    for k in ("abc_0.4_2.0", "abc_0.4_1.0", "tail_rejected_F_le_10", "abc_and_tail"):
        rows.append([v, k, d[k]["frac"], d[k]["P1_success"], d[k]["P1_median"]])
w("S8b_constrained_restrictions", ["variant", "restriction", "fraction_retained", "P1_success", "P1_median_F"], rows)
D = J("decomp"); rows = []
for m, d in D.items():
    for k, s in d.items():
        rows.append([m, k, s["median"], s["p95"], s["gain"]])
w("S9_matrix_only_decomposition", ["model", "arm", "median_F", "p95_F", "gain_over_peptide_alone"], rows)
R = J("revers")["ctox_ratio"]; rows = []
for mult, v in R.items():
    for arm in ("R1_reference_PE", "P1_dual_lock", "P5_no_self_limiting_PE"):
        a = v[arm]
        for key, val in a["grid"].items():
            dh, tr = key.split("|"); rows.append([v["C_tox_over_EC50_median"], arm, dh, tr, val])
w("S10_reversibility_grid", ["Ctox_over_EC50_median", "arm", "DH_max_h", "recovery_max_h", "P_success"], rows)
O = J("optim"); rows = []
for m, v in O.items():
    for k in ("nominal_train", "nominal_test", "optimized_train", "optimized_test", "reference_test"):
        if k in v:
            s = v[k]; rows.append([m, k, json.dumps(v.get("best", {})), s["p_success"], s["p_safe"], s["median"], s["p95"], s["frac_gt_10"]])
w("S11_optimization", ["model", "set", "design_best", "P_success", "P_safe", "median_F", "p95_F", "frac_F_gt_10"], rows)
