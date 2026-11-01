# -*- coding: utf-8 -*-
"""Sections 1-4 of the manuscript (design-dependent text; result numbers are inserted by build_manuscript.py in sections 5-9). Block kinds as in build_infra."""
import oral_peptide_mc as M

B = []
add = B.append
PR = M.PRIORS
A = M.ARMS


def sig(x):
    return ("%.2f" % x).rstrip("0").rstrip(".")


add(("h1", "1. Introduction"))
add(("p", "Peptide drugs with molecular masses of 2 to 4 kDa combine high potency with good selectivity, but almost all are still injected. After oral administration they meet three barriers in series: "
          "luminal and brush-border proteolysis, low passive permeability across the epithelium, and, for some sequences, apical efflux [@Maher2016;@Drucker2020]. Formulations built on "
          "permeation enhancers (PEs) are the most advanced response. They have delivered one approved tablet and several late-stage candidates, but absolute bioavailability has stayed at about 1% "
          "and is highly variable [@Brayden2020]. For oral semaglutide co-formulated with the enhancer salcaprozate sodium (SNAC), a population pharmacokinetic analysis gave a bioavailability of "
          "0.8% under the recommended dosing conditions, with a within-subject variability of 137% [@Overgaard2021]. A higher and more reproducible bioavailability would lower dose, cost and pill burden."))
add(("p", "The mechanistic work on SNAC points to an unusual feature of enhancer tablets: the effect is local. Absorption of semaglutide was transcellular, confined to the area next to the tablet "
          "surface, only transiently enhanced, and dependent on local buffering that protected the peptide from enzymatic degradation [@Buckley2018]. A comparison of SNAC with sodium caprate "
          "showed that both give single-digit, variable increases in bioavailability and that the mechanisms at in vivo doses remain debated [@Twarog2019]. If efficacy is a property of a micro-niche in which "
          "peptide, enhancer and protection coincide in space and time, then the controllable quantities are the size of the niche, how long it persists, and how closely the release of its "
          "components is synchronized. Conventional enhancer tablets do not control any of these deliberately, and no quantitative framework links them to bioavailability and to the reversibility of the epithelial effect."))
add(("p", "Efflux is the second cell-level barrier. D-$\\alpha$-tocopheryl polyethylene glycol 1000 succinate (TPGS) inhibits P-glycoprotein in Caco-2 monolayers without being a substrate itself [@Collnot2006] "
          "and has been reported to raise the oral bioavailability of P-glycoprotein substrates [@Dintaman1999], whereas piperine inhibits P-glycoprotein and CYP3A4 and carries a drug-interaction liability [@Bhardwaj2002]. "
          "Whether a 2 to 4 kDa hydrophilic peptide is an efflux substrate is peptide-specific, so the value of a modulator is conditional. Hydrophobic ion pairing is a non-covalent route to higher "
          "lipophilicity that is retained here only as an optional enabler [@Griesser2017;@Ristroph2019]."))
add(("p", "Here we ask a design question that experiments cannot yet answer: under which conditions would a co-localized, reversible enhancer and efflux-modulator architecture, which we call Dual-Lock, "
          "raise oral bioavailability above 4% without an unacceptable epithelial exposure, and under which conditions would it fail? We answer it with a mechanistic Monte Carlo analysis whose inputs are "
          "evidence-graded priors, and we treat the answer as a set of falsifiable predictions rather than a performance claim. The study makes five contributions."))
for t in ["a mechanistic model with every equation, unit, initial condition and numerical method stated (Section 3);",
          "a paired Monte Carlo design in which all arms share the same random draws, so that ablations and negative controls isolate co-localization, enhancer clearance, protection and efflux modulation (Sections 4 and 5);",
          "global variance-based sensitivity analysis and a decomposition of output uncertainty into peptide, enhancer, modulator and formulation contributions (Section 5);",
          "a benefit-risk (Pareto) analysis and falsification-region maps that state, in the units of future experiments, where the architecture fails (Sections 5 and 6);",
          "a structural-uncertainty analysis that replaces the lumped niche with a spatial diffusion-reaction model, and a comparison of the model with published human observations that were not used to calibrate it (Section 6)."]:
    add(("bullet", t))
add(("p", "The work is entirely computational. It contains no experimental data, none of its probabilities is an empirical probability of success, and the reference arm is calibrated to a published baseline, so the "
          "reference result is not a validation. Section 8 gives the experiments that would test the predictions."))

add(("h1", "2. Mechanistic framework"))
add(("p", "Dual-Lock is an excipient-only co-formulation of the unmodified peptide, delivered as a conventional solid dosage form, that combines four modules."))
add(("bullet", "**M1, co-localization matrix.** A mucoadhesive, pH-triggered matrix releases peptide, enhancer and modulator in step into a thin, slowly diluted niche against the mucosa. It is described by the niche thickness $\\delta$, the wash-out rate $k_{dil}$, the release rate $k_{rel}$ and the release mismatch $m$."))
add(("bullet", "**M2, reversible, self-limiting enhancer.** A transcellular amphiphile is selected for fast partition-and-departure kinetics and enzymatic or pH-driven clearance ($k_{dep}$), and is dosed at half the reference amount. Its effect on permeability follows a Hill relation in the local concentration."))
add(("bullet", "**M3, efflux modulator.** A pharmaceutical-grade P-glycoprotein modulator (TPGS first-line) lowers the apical efflux ratio, but only for peptides that are efflux substrates."))
add(("bullet", "**M4, local buffer and protease shield.** A buffering and protease-protection module lowers the effective proteolysis rate by the protection factor $Pr$."))
add(("p", "The framework gives four mechanistic hypotheses that the model is designed to test: **H1**, confinement of the niche raises the absorbed fraction of peptide at a lower total enhancer dose; **H2**, a self-limiting enhancer is the price of "
          "reversibility, because confinement also raises the local enhancer concentration; **H3**, efflux modulation matters only for efflux substrates; and **H4**, protection adds a smaller benefit than confinement."))
add(("p", "*What the niche thickness is.* In the model, $\\delta$ is the volume-to-contact-area ratio of the confined compartment, $\\delta=V/A_m$ (cm). It is a design target for the matrix, not a measured mucus thickness. "
          "The adherent mucus gel measured in rat intestine is 20 to 830 $\\mu$m thick (jejunum about 123 $\\mu$m, colon about 830 $\\mu$m) [@Atuma2001], that is, 0.002 to 0.083 cm, whereas the nominal Dual-Lock value of 0.30 cm "
          "is the thickness of a continuously mixed fluid-gel layer holding the dose, and the reference value of 1.0 cm represents the larger volume around a conventional tablet. Because transport across a layer of this "
          "thickness by molecular diffusion alone would take hours (Section 6.1), the lumped model implicitly assumes mixing within the niche; the spatial model of Section 3.5 tests that assumption."))

add(("h1", "3. Mathematical methods"))
add(("h2", "3.1 State variables and balances"))
add(("p", "Quantities are expressed per dose. The peptide dose is 1 (all peptide amounts are fractions of the dose); the enhancer and modulator doses $E_0$ and $I_0$ are in mmol. The peptide pools are the tablet pool $P_t$, the amount in the niche $P$, "
          "and the cumulative amounts absorbed ($Q_{abs}$), degraded ($Q_{deg}$) and washed out ($Q_{dil}$). Release from the tablet is first order, with the enhancer and modulator released at the same rate scaled by the mismatch factor $m$:"))
add(("eq", r"\frac{dP_t}{dt}=-k_{rel}P_t,\qquad \frac{dE_t}{dt}=-m\,k_{rel}E_t,\qquad \frac{dI_t}{dt}=-m\,k_{rel}I_t", 1))
add(("p", "The niche is a single compartment of volume $V=A_m\\delta$ (mL, with $A_m$ in cm$^2$ and $\\delta$ in cm). With the apical uptake rate constant $k_{abs}$ defined below,"))
add(("eq", r"\frac{dP}{dt}=k_{rel}P_t-\left[k_{dil}+k_{deg}(1-Pr)+k_{abs}\right]P", 2))
add(("eq", r"\frac{dE}{dt}=m\,k_{rel}E_t-\left(k_{dil}+k_{dep}\right)E,\qquad \frac{dI}{dt}=m\,k_{rel}I_t-k_{dil}\,I", 3))
add(("p", "All state variables start at zero except $P_t=1$, $E_t=E_0$ and $I_t=I_0$. The cumulative amounts obey $dQ_{abs}/dt=k_{abs}P$, $dQ_{deg}/dt=k_{deg}(1-Pr)P$ and $dQ_{dil}/dt=k_{dil}P$, so that the peptide balance"))
add(("eq", r"P_t+P+Q_{abs}+Q_{deg}+Q_{dil}=1", 4))
add(("p", "holds at all times and serves as a numerical check (maximum error below $10^{-9}$ in all arms)."))
add(("h2", "3.2 Enhancement, efflux, absorption and exposure"))
add(("p", "The enhancer and modulator concentrations in the niche, in mM, are $C_E=10^3E/V$ and $C_I=10^3I/V$. Enhancement of transcellular permeability follows a Hill relation,"))
add(("eq", r"H=\frac{C_E^{\,h}}{EC_{50}^{\,h}+C_E^{\,h}},\qquad P_{eff}=P_0\left(1+E_{max}H\right)", 5))
add(("p", "where $P_0$ is the baseline apparent permeability of the peptide (cm h$^{-1}$), $EC_{50}$ the enhancer concentration for half-maximal enhancement, $E_{max}$ the maximal fold-increase and $h$ the Hill coefficient. "
          "For a peptide that is an efflux substrate the modulator lowers the efflux ratio from its baseline value $ER_0$ towards 1; for a non-substrate $ER=1$:"))
add(("eq", r"ER=1+\frac{ER_0-1}{1+C_I/IC_{50}},\qquad k_{abs}=\frac{P_{eff}}{\delta\,ER}", 6))
add(("p", "The absorbed fraction and the bioavailability follow from the post-absorptive survival $S_{post}$ (intra-enterocyte handling and first pass),"))
add(("eq", r"F=S_{post}\,Q_{abs}(T),\qquad Q_{abs}(T)=\int_0^{T}k_{abs}P\,dt,\qquad T=8\ \mathrm{h}", 7))
add(("p", "Epithelial exposure to the enhancer is summarized by damage-hours, using a steep Hill function of the enhancer concentration with half-damage concentration $C_{tox}$,"))
add(("eq", r"D(t)=\frac{C_E^{\,3}}{C_{tox}^{\,3}+C_E^{\,3}},\qquad DH=\int_0^{T}D(t)\,dt", 8))
add(("p", "Reversibility requires both a bounded exposure and a bounded duration. With $t_{clear}$ the last time at which $C_E>0.1\\,EC_{50}$ and $t_{rev}$ the epithelial recovery time after the enhancer has fallen below that threshold,"))
add(("eq", r"\mathrm{reversible}\iff DH\le 0.5\ \mathrm{h}\ \ \mathrm{and}\ \ t_{clear}+t_{rev}\le 4\ \mathrm{h};\qquad \mathrm{success}\iff F\ge 0.04\ \mathrm{and\ reversible}", 9))
add(("p", "The thresholds of 4% bioavailability and reversibility within hours are those of the published challenge that motivated the work; they are design targets, not regulatory limits. Success is therefore a model-conditional indicator, "
          "and its mean over the prior is called P(success)."))
add(("p", "In this paper \"reversible\" always means *predicted exposure reversibility*, a model-defined indicator, and never demonstrated safety. Its two parts correspond to measurable epithelial endpoints: the time $t_{clear}+t_{rev}$ stands for the recovery of transepithelial electrical resistance (TEER) after the enhancer has fallen below its effective concentration, "
          "and damage-hours stand for cumulative exposure to concentrations that reduce viability (for example lactate dehydrogenase release). Both are surrogates with known limitations. In isolated rat mucosae, SNAC and sodium caprate lowered TEER and caused minor epithelial damage at the concentrations that increased fluxes [@Twarog2021], "
          "so the model's assumption that damage occurs only at about three times the effective concentration is tested explicitly (Section 5.9)."))  # REV2
add(("h2", "3.3 Calibration"))
add(("p", "Baseline permeability is the one parameter fixed by data. The median of $P_0$ is set by bisection so that the median bioavailability of the reference arm R1 (a SNAC-like enhancer tablet with a 1.0 cm niche) equals 1.2%, the mid-point of the 1 to 2% range usually quoted for enhancer tablets "
          "[@Brayden2020]:"))
add(("eq", r"\mathrm{median}\left\{F_{R1}(P_0)\right\}=0.012", 10))
add(("p", "The calibrated median is $P_0=0.00975$ cm h$^{-1}$ (about $2.7\\times10^{-6}$ cm s$^{-1}$). Because the reference arm is fitted, its bioavailability is not a validation; only the contrasts between arms, which share all other priors, carry information."))
add(("h2", "3.4 Numerical method and verification"))
add(("p", "Equations (1) to (3) and the cumulative terms are integrated with the explicit Euler method at $\\Delta t=0.01$ h over 8 h, vectorized over draws. Halving the time step four times changed the median bioavailability by less than 2% (Section 5.3 and Supplementary Fig. S1), "
          "and extending the horizon to 24 h changed it by less than the stated tolerance (Section 5.3). Verification consists of the peptide balance of Eq. (4), module-logic tests (the efflux module cannot lower $F$ and has no effect on non-substrates; the enhancer raises $F$ over peptide alone), "
          "determinism for a fixed seed, and the tests of the extension listed in Supplementary Table S6. An early prototype contained a units error in the concentration (a factor of $10^3$) that made the enhancer and modulator appear inert; it was found by the module-logic tests, and all results here use the corrected code. "
          "These tests verify the software, not the biology."))
add(("h2", "3.5 Alternative structural models"))
add(("p", "The lumped niche assumes instantaneous mixing. Two alternatives test that assumption. In the *stagnant-layer model* a serial diffusion resistance is added between the mixed niche and the epithelial membrane, replacing $P_{eff}$ in Eq. (6) by"))
add(("eq", r"P_{eff,s}=\left(\frac{1}{P_{eff}}+\frac{h_s}{D_m}\right)^{-1}", 11))
add(("p", "with a layer thickness $h_s$ and a peptide diffusivity $D_m$ in the layer. In the *spatial model* the niche is resolved along the coordinate $x$ from the matrix and lumen side ($x=0$) to the epithelial surface ($x=\\delta$). The peptide concentration $c_P(x,t)$ (fraction of dose per mL) obeys"))
add(("eq", r"\frac{\partial c_P}{\partial t}=D_P\frac{\partial^2 c_P}{\partial x^2}-k_{deg}(1-Pr)\,c_P", 12))
add(("p", "with a wash-out flux at $x=0$ through the mass-transfer coefficient $h_0=k_{dil}\\delta$ and the release flux $k_{rel}P_t/A_m$ entering there, and an absorbing boundary at $x=\\delta$ with the same permeability as the lumped model:"))
add(("eq", r"D_P\left.\frac{\partial c_P}{\partial x}\right|_{0}=h_0\,c_P(0,t)-\frac{k_{rel}P_t}{A_m},\qquad -D_P\left.\frac{\partial c_P}{\partial x}\right|_{\delta}=\frac{P_{eff}}{ER}\,c_P(\delta,t)", 13))
add(("p", "The enhancer and modulator obey the same equation with their own diffusivities, with the clearance $k_{dep}$ acting on the enhancer only; $P_{eff}$ and $ER$ are evaluated from the enhancer and modulator concentrations in the epithelial cell. "
          "The scheme is a backward-Euler finite-volume discretization (24 cells, $\\Delta t=0.02$ h), unconditionally stable and exactly conservative, with 1 to $10^6$ times the molecular diffusivity (an effective-dispersion multiplier that stands for luminal mixing). "
          "In the limit of strong mixing it reduces to the lumped model (tested; Section 5.3), and mesh and time-step refinement changed the median bioavailability by less than 3%. For each structural variant the baseline permeability is recalibrated so that R1 again has a median of 1.2%, "
          "because the calibration anchor must hold in every structure."))

add(("p", "*Uptake-constrained model.* The compartment model lets the uptake rate constant $P_{eff}/(\\delta\\,ER)$ grow without bound as the niche is thinned, which is the origin of the large matrix-only gain discussed in Section 5.8. Two physical limits are therefore added. "
          "The stagnant layer of Eq. (11) caps the permeability that the epithelium can express. The dissolved peptide in the niche cannot exceed its solubility $S$ (mg mL$^{-1}$); with the dose $D$ (mg) and the niche volume $V$,"))
add(("eq", r"P\le P_{max}=\frac{S\,V}{D},\qquad \frac{dP_t}{dt}=-k_{rel}P_t\left(1-\frac{P}{P_{max}}\right)_{+}", 14))
add(("p", "so that the peptide released into a small niche can only dissolve up to its solubility and the remainder stays in the tablet pool (the balance of Eq. (4) is unchanged). The dose and the solubility are grade-D priors (Table 2). "
          "Both limits are monotone: they can only lower the bioavailability of a given draw (tested), and the baseline permeability is recalibrated for every variant so that R1 keeps a median of 1.2%. Peptide-enhancer complexation is not modelled; measured affinities for exenatide with SNAC and sodium caprate are low (dissociation constants of about 10 to 100 $\\mu$M) [@Twarog2022]."))
add(("h1", "4. Monte Carlo methods"))
add(("h2", "4.1 Priors, evidence grades and correlations"))
add(("p", "Each uncertain input has a prior and an evidence grade: A, a directly applicable measurement; B, a measurement for a close analogue; C, an order-of-magnitude literature range; D, an author assumption without direct source. No input reaches grade A. "
          "Table 2 lists the common priors (log-normal unless stated, with the median and the standard deviation of the natural logarithm, $\\sigma_{\\ln}$). The probability that the peptide is an efflux substrate is a scenario axis, not a prior, "
          "and is set to 0.5 in the base case."))
pri = []
rows_pr = [("$P_0$", "baseline apparent permeability", "cm h$^{-1}$", "C", "0.00975 (calibrated), $\\sigma_{\\ln}$ 0.60"),
           ("$k_{deg}$", "luminal proteolysis rate without protection", "h$^{-1}$", "D", "2.0, 0.50"),
           ("$S_{post}$", "survival after apical entry", "-", "D", "beta, mean 0.55, sd 0.15"),
           ("$ER_0$", "efflux ratio if the peptide is a substrate", "-", "C", "2.5, 0.50"),
           ("$A_m$", "mucosal contact area of the niche", "cm$^2$", "D", "30, 0.30"),
           ("$EC_{50}$", "enhancer concentration for half-maximal effect", "mM", "C", "25, 0.50"),
           ("$E_{max}$", "maximal fold-increase of permeability", "-", "D", "80, 0.50"),
           ("$h$", "Hill coefficient of the enhancer effect", "-", "D", "2.0, 0.25"),
           ("$C_{tox}$", "enhancer concentration for half-maximal damage", "mM", "D", "80, 0.40"),
           ("$t_{rev}$", "epithelial recovery time", "h", "D", "1.0, 0.50"),
           ("$IC_{50}$", "efflux-modulator potency", "mM", "C", "0.5, 0.80"),
           ("$D$", "peptide dose (uptake-constrained model)", "mg", "D", "100, 0.35"),
           ("$S$", "peptide solubility (uptake-constrained model)", "mg mL$^{-1}$", "D", "20, 0.80"),
           ("$h_s$", "stagnant-layer thickness", "cm", "C", "0.015, 0.70 (20 to 830 $\\mu$m)"),
           ("$D_m$", "peptide diffusivity in the layer", "cm$^2$ s$^{-1}$", "C", "$10^{-6}$, 0.50")]
add(("table", "**Table 2.** Common priors (log-normal: median, $\\sigma_{\\ln}$) with evidence grade. Bases: $P_0$ calibrated to the 1 to 2% enhancer-tablet baseline [@Brayden2020]; $ER_0$ from Caco-2 efflux ratios; $EC_{50}$ in the millimolar range reported for SNAC and sodium caprate [@Maher2016;@Twarog2019]; "
          "$IC_{50}$ spans piperine at 15 to 74 $\\mu$M in Caco-2 [@Bhardwaj2002] and sub-millimolar TPGS [@Collnot2006].",
     ["Symbol", "Meaning", "Unit", "Grade", "Prior"], [list(r) for r in rows_pr], [0.8, 3.0, 0.9, 0.6, 2.3], 16))
add(("p", "Formulation-controlled inputs are log-normal around arm-specific medians (dispersion $\\sigma_{\\ln}$ 0.30 for $\\delta$, $k_{dil}$ and $k_{dep}$, 0.20 for $k_{rel}$, and the arm's mismatch value for $m$) and beta-distributed for $Pr$ (sd 0.10). "
          "Table 3 lists the arms. **All arm settings of the matrix are design hypotheses of grade D.** The Dual-Lock enhancer has the same $EC_{50}$, $E_{max}$ and toxicity priors as the reference enhancer, so that no benefit can arise from assumed potency."))
arm_rows = []
nm = [("R0_peptide_alone", "R0 peptide alone"), ("R1_reference_PE", "R1 reference PE (SNAC-like)"), ("R2_reference_PE_efflux", "R2 R1 + efflux module"), ("P1_dual_lock", "P1 Dual-Lock"),
      ("P2_no_efflux_module", "P2 P1 without modulator"), ("P3_no_colocalisation", "P3 P1 without co-localization"), ("P4_no_protection", "P4 P1 with reference protection"),
      ("P5_no_self_limiting_PE", "P5 P1 without self-limiting PE")]
for k, lab in nm:
    a = A[k]; arm_rows.append([lab, "%g" % a["pe"], "%g" % a["inh"], "%g" % a["delta"], "%g" % a["k_dil"], "%g" % a["mism"], "%g" % a["prot"], "%g" % a["k_dep"], "%g" % a["k_rel"]])
add(("table", "**Table 3.** Arm definitions (medians of formulation-controlled inputs; all grade D). PE and modulator doses in mmol; $\\delta$ in cm; $k_{dil}$, $k_{dep}$, $k_{rel}$ in h$^{-1}$; $m$ is the log-scale dispersion of the release mismatch; $Pr$ the mean protection factor.",
     ["Arm", "PE", "Modulator", "$\\delta$", "$k_{dil}$", "$m$", "$Pr$", "$k_{dep}$", "$k_{rel}$"], arm_rows, [3.4, 0.7, 1.0, 0.7, 0.8, 0.6, 0.6, 0.8, 0.8], 15))
add(("p", "In the base case the inputs are independent. Biological and pharmacological parameters are in fact coupled, so correlated priors are added through a Gaussian copula with correlation matrix $R$ on the normal scores: "
          "$U=\\Phi(Z)$, $Z\\sim N(0,R)$, followed by the inverse marginal distributions [@Iman1982]. The assumed (grade D) correlations are $\\rho(EC_{50},C_{tox})=+0.6$ (potent enhancers also damage the epithelium at lower concentration), "
          "$\\rho(EC_{50},E_{max})=-0.3$ and $\\rho(C_{tox},t_{rev})=-0.4$ (more damaging agents recover more slowly). The central coupling is varied from 0 to 0.9 as a sensitivity analysis."))
add(("h2", "4.2 Sampling, paired design and convergence"))
add(("p", "Each arm is run with 20,000 draws and a fixed seed. For contrasts between arms all arms use the same matrix of uniform variates, that is, common random numbers, so that biological parameters are identical across arms within a draw and paired differences remove the "
          "Monte Carlo noise of the biology. The sampler was re-implemented to map uniform variates to the priors; it reproduced the success probability of the original independent sampler within sampling error (Section 5.1). "
          "Convergence was assessed from running estimates up to $10^5$ draws, batch means, 20 independent replicate runs, and a time-step and horizon study."))
add(("h2", "4.3 Global sensitivity and uncertainty decomposition"))
add(("p", "Global uncertainty and sensitivity analysis follows established practice in systems biology [@Marino2008]. Variance-based sensitivity indices were estimated for $\\log_{10}F$, the success indicator and $\\log_{10}DH$ with the Saltelli design (a scrambled Sobol sequence, $N=2^{13}$ base samples and $N(d+2+g)$ model evaluations for $d=18$ inputs and $g=4$ groups). "
          "The first-order index is $S_i=V[E(Y|X_i)]/V(Y)$ and the total-order index $S_{Ti}=E[V(Y|X_{\\sim i})]/V(Y)$, estimated with the Saltelli and Jansen estimators [@Saltelli2010;@Jansen1999;@Sobol2001]; intervals are 95% bootstrap intervals (200 resamples). "
          "Group indices combine the inputs into peptide and epithelium ($P_0$, $k_{deg}$, $S_{post}$, $ER_0$, substrate status, $A_m$), enhancer pharmacology ($EC_{50}$, $E_{max}$, $h$, $C_{tox}$, $t_{rev}$), the modulator ($IC_{50}$) and the formulation architecture ($\\delta$, $k_{dil}$, $m$, $Pr$, $k_{dep}$, $k_{rel}$). "
          "The estimator was verified on the Ishigami function (Supplementary Table S6)."))
add(("h2", "4.4 Ablations, negative controls, benefit-risk and falsification analyses"))
add(("p", "Ablations remove or add one component at a time (removal from P1; addition to R1), with the paired difference in P(success) as the effect measure. Negative controls are the matrix without enhancer or modulator (M0), the reference enhancer at half dose, and "
          "the architecture at full dose. The benefit-risk analysis sweeps the enhancer dose from 0.05 to 3 mmol for the reference, Dual-Lock, and the two ablations that matter for reversibility, and plots P($F\\ge4\\%$) against the probability that the exposure is not reversible. "
          "The falsification map evaluates P(success) on a grid of niche thickness and enhancer dose (and of niche thickness and wash-out), with a predefined benefit region B: median $F\\ge4\\%$ and P(reversible) $\\ge0.80$ for some dose."))
add(("h2", "4.5 Calibration-independent checks"))
add(("p", "Because the reference arm is fitted, three checks do not rely on the calibration. First, an observation-conditioned analysis (approximate Bayesian computation by rejection [@Beaumont2002]) keeps only the draws in which R1 falls in a published range "
          "(0.4 to 2.0%, or the 0.4 to 1.0% of oral semaglutide) and reports the Dual-Lock outcome among them. Second, the spread, the upper tail and the qualitative dependences of the reference arm are compared with published human observations "
          "[@Overgaard2021;@Granhall2019;@Brayden2020;@Pechenov2021;@Buckley2018]. Third, the structural variants of Section 3.5 are recalibrated and compared. These checks can reveal inconsistency; they cannot validate the Dual-Lock arm, which no data address."))
add(("h2", "4.6 Software and reproducibility"))
add(("p", "The model (oral_peptide_mc.py, unchanged from the first version apart from tie-corrected rank correlation), the extension (dualock_lib.py, run_extended.py, make_figures.py), the automated tests, all numerical results (JSON), the parameter provenance table, and the figure scripts are in a versioned repository "
          "archived with a DOI (Data availability). All numbers in the text are generated from the result files."))
add(("h2", "4.7 Analyses added in the revision"))
add(("p", "*Uptake constraints and the matrix-only control.* The reference arm, the matrix-only control (P1 architecture without enhancer or modulator), peptide alone, Dual-Lock and its main ablations were run in five variants: the original lumped model; the stagnant layer; the solubility limit; both limits (the uptake-constrained model); and the uptake-constrained model in which damage occurs at the efficacious concentration ($C_{tox}$ set to the median $EC_{50}$). "
          "Each variant was recalibrated to the same 1.2% anchor. Within each variant the Dual-Lock outcome was also evaluated among the draws in which R1 lies in the published ranges, among the draws in which no arm exceeds 10% (a bounding analysis that treats larger values as non-credible), and under both restrictions."))
add(("p", "*Reversibility endpoints.* The success indicator was recomputed on a grid of damage-hours limits (0.5 to 0.001 h) and recovery limits (2 to 6 h), and with the median $C_{tox}$ reduced to 30%, 10% and 3% of its base value (potency-to-toxicity ratio from 3.2 down to 0.1)."))
add(("p", "*Optimized design out of sample.* The release rate (1.2 to 4 h$^{-1}$), enhancer dose (0.35 to 1.0 mmol) and wash-out (0.4 to 1.5 h$^{-1}$) of the nominal design were varied on a grid of 72 combinations; the combination with the highest P(success) on 8,000 training draws defines the optimized design, which, together with the nominal design and the reference, was evaluated on 20,000 independent test draws (new seed). "
          "The difference between the training and the test estimate is the optimism of the optimization."))
