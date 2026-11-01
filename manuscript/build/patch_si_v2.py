# -*- coding: utf-8 -*-
import os
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_si.py")
s = open(p, encoding="utf-8").read()
if "REV2" in s:
    raise SystemExit("already")
s = s.replace("16 tests). They verify", "19 tests). They verify").replace("`python run_extended.py --test`: 16 tests", "python run_extended.py --test: 19 tests")
add_tests = '''         ("Solubility limit and input modifiers", "the solubility limit never raises $F$ for the same draws and the mass balance closes to $10^{-9}$; lowering $C_{tox}$ never lowers damage-hours."),
'''
s = s.replace('         ("Determinism and helpers"', add_tests + '         ("Determinism and helpers"', 1)
new = '''
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
'''
s = s.replace('add(("refs",))', new + '\nadd(("refs",))', 1)
open(p, "w", encoding="utf-8").write(s)
print("patched")
