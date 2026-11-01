# -*- coding: utf-8 -*-
"""Revision patch of text_front.py (methods additions for the uptake-constrained model, reversibility endpoints and out-of-sample optimization)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
p = os.path.join(HERE, "text_front.py")
s = open(p, encoding="utf-8").read()
if "REV2" in s:
    raise SystemExit("already patched")

# (a) reversibility wording after Eq. (9)
a = 'add(("h2", "3.3 Calibration"))'
b = '''add(("p", "In this paper \\"reversible\\" always means *predicted exposure reversibility*, a model-defined indicator, and never demonstrated safety. Its two parts correspond to measurable epithelial endpoints: the time $t_{clear}+t_{rev}$ stands for the recovery of transepithelial electrical resistance (TEER) after the enhancer has fallen below its effective concentration, "
          "and damage-hours stand for cumulative exposure to concentrations that reduce viability (for example lactate dehydrogenase release). Both are surrogates with known limitations. In isolated rat mucosae, SNAC and sodium caprate lowered TEER and caused minor epithelial damage at the concentrations that increased fluxes [@Twarog2021], "
          "so the model's assumption that damage occurs only at about three times the effective concentration is tested explicitly (Section 5.9)."))  # REV2
add(("h2", "3.3 Calibration"))'''
assert a in s; s = s.replace(a, b, 1)

# (b) uptake-constrained model after the spatial paragraph (before section 4)
a = 'add(("h1", "4. Monte Carlo methods"))'
b = '''add(("p", "*Uptake-constrained model.* The compartment model lets the uptake rate constant $P_{eff}/(\\\\delta\\\\,ER)$ grow without bound as the niche is thinned, which is the origin of the large matrix-only gain discussed in Section 5.8. Two physical limits are therefore added. "
          "The stagnant layer of Eq. (11) caps the permeability that the epithelium can express. The dissolved peptide in the niche cannot exceed its solubility $S$ (mg mL$^{-1}$); with the dose $D$ (mg) and the niche volume $V$,"))
add(("eq", r"P\\le P_{max}=\\frac{S\\,V}{D},\\qquad \\frac{dP_t}{dt}=-k_{rel}P_t\\left(1-\\frac{P}{P_{max}}\\right)_{+}", 14))
add(("p", "so that the peptide released into a small niche can only dissolve up to its solubility and the remainder stays in the tablet pool (the balance of Eq. (4) is unchanged). The dose and the solubility are grade-D priors (Table 2). "
          "Both limits are monotone: they can only lower the bioavailability of a given draw (tested), and the baseline permeability is recalibrated for every variant so that R1 keeps a median of 1.2%. Peptide-enhancer complexation is not modelled; measured affinities for exenatide with SNAC and sodium caprate are low (dissociation constants of about 10 to 100 $\\\\mu$M) [@Twarog2022]."))
add(("h1", "4. Monte Carlo methods"))'''
assert a in s; s = s.replace(a, b, 1)

# (c) Table 2 rows
a = '           ("$IC_{50}$", "efflux-modulator potency", "mM", "C", "0.5, 0.80")]'
b = '''           ("$IC_{50}$", "efflux-modulator potency", "mM", "C", "0.5, 0.80"),
           ("$D$", "peptide dose (uptake-constrained model)", "mg", "D", "100, 0.35"),
           ("$S$", "peptide solubility (uptake-constrained model)", "mg mL$^{-1}$", "D", "20, 0.80"),
           ("$h_s$", "stagnant-layer thickness", "cm", "C", "0.015, 0.70 (20 to 830 $\\\\mu$m)"),
           ("$D_m$", "peptide diffusivity in the layer", "cm$^2$ s$^{-1}$", "C", "$10^{-6}$, 0.50")]'''
assert a in s; s = s.replace(a, b, 1)

# (d) new methods subsection at the end of section 4
s = s.rstrip() + '''
add(("h2", "4.7 Analyses added in the revision"))
add(("p", "*Uptake constraints and the matrix-only control.* The reference arm, the matrix-only control (P1 architecture without enhancer or modulator), peptide alone, Dual-Lock and its main ablations were run in five variants: the original lumped model; the stagnant layer; the solubility limit; both limits (the uptake-constrained model); and the uptake-constrained model in which damage occurs at the efficacious concentration ($C_{tox}$ set to the median $EC_{50}$). "
          "Each variant was recalibrated to the same 1.2% anchor. Within each variant the Dual-Lock outcome was also evaluated among the draws in which R1 lies in the published ranges, among the draws in which no arm exceeds 10% (a bounding analysis that treats larger values as non-credible), and under both restrictions."))
add(("p", "*Reversibility endpoints.* The success indicator was recomputed on a grid of damage-hours limits (0.5 to 0.001 h) and recovery limits (2 to 6 h), and with the median $C_{tox}$ reduced to 30%, 10% and 3% of its base value (potency-to-toxicity ratio from 3.2 down to 0.1)."))
add(("p", "*Optimized design out of sample.* The release rate (1.2 to 4 h$^{-1}$), enhancer dose (0.35 to 1.0 mmol) and wash-out (0.4 to 1.5 h$^{-1}$) of the nominal design were varied on a grid of 72 combinations; the combination with the highest P(success) on 8,000 training draws defines the optimized design, which, together with the nominal design and the reference, was evaluated on 20,000 independent test draws (new seed). "
          "The difference between the training and the test estimate is the optimism of the optimization."))
'''
open(p, "w", encoding="utf-8").write(s)
print("patched")
