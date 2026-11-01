# -*- coding: utf-8 -*-
"""make_figures.py -- figures for the Dual-Lock manuscript from results_ext/*.json (run after run_extended.py)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "results_ext"
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130, "savefig.dpi": 300, "font.family": "DejaVu Sans"})
J = lambda n: json.load(open(os.path.join(OUT, "ext_%s.json" % n)))
C = {"R": "#6b6b6b", "P": "#1f77b4", "Q": "#d62728", "G": "#2ca02c", "O": "#ff7f0e"}
LAB = {"R0_peptide_alone": "R0 peptide alone", "R1_reference_PE": "R1 reference PE", "R2_reference_PE_efflux": "R2 ref. PE + efflux", "P1_dual_lock": "P1 Dual-Lock",
       "P2_no_efflux_module": "P2 no efflux module", "P3_no_colocalisation": "P3 no co-localization", "P4_no_protection": "P4 ref. protection",
       "P5_no_self_limiting_PE": "P5 no self-limiting PE", "M0_matrix_only": "M0 matrix only (no PE, no modulator)", "R1_half_dose": "R1 at half PE dose (0.5 mmol)",
       "P1_full_dose": "P1 at full PE dose (1.0 mmol)", "P1_minus_delta": "P1, niche thickness reverted", "P1_minus_k_dil": "P1, wash-out reverted",
       "P1_minus_mism": "P1, release mismatch reverted", "P1_minus_k_rel": "P1, release rate reverted", "R1_plus_confinement_bundle": "R1 + confinement bundle",
       "R1_plus_delta": "R1 + thin niche", "R1_plus_k_dil": "R1 + slow wash-out", "R1_plus_mism": "R1 + synchronous release", "R1_plus_k_rel": "R1 + slower release",
       "R1_plus_prot": "R1 + protection 0.70", "R1_plus_k_dep": "R1 + fast PE clearance", "R1_plus_halfPE_efflux": "R1 + half PE + efflux"}


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), bbox_inches="tight"); plt.close(fig); print("wrote", name)


def fig_schematic():
    fig, ax = plt.subplots(figsize=(7.2, 3.9)); ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")
    def box(x, y, w, h, fc, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=fc, ec=ec, lw=1.0))
    def arr(x1, y1, x2, y2, c="k"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10, lw=1.1, color=c, shrinkA=0, shrinkB=0))
    box(1, 30, 19, 20, "#fff4e0", "#b36b00")
    ax.text(10.5, 46.5, "Tablet pools", ha="center", fontsize=7.5, weight="bold")
    for k, t in enumerate(["$P_t$ peptide", "$E_t$ enhancer", "$I_t$ modulator"]):
        ax.text(10.5, 41.5 - 4.2 * k, t, ha="center", fontsize=7.5)
    arr(20.5, 40, 31.5, 40); ax.text(26, 43.2, "release $k_{rel}$", ha="center", fontsize=6.8); ax.text(26, 36.3, "mismatch $m$", ha="center", fontsize=6.8)
    box(32, 22, 38, 30, "#e8f4ea", "#2e7d32")
    ax.text(51, 48.3, "Co-localized niche, $V=A_m\delta$", ha="center", fontsize=8, weight="bold")
    ax.text(51, 43.0, "$P$ peptide    $E$ enhancer    $I$ modulator", ha="center", fontsize=7.2)
    ax.text(51, 37.5, "Hill enhancement $H(C_E)$", ha="center", fontsize=7.2)
    ax.text(51, 33.0, "efflux ratio $ER(C_I)$", ha="center", fontsize=7.2)
    ax.text(51, 27.3, "protection $Pr$:  proteolysis $k_{deg}(1-Pr)$", ha="center", fontsize=7.2)
    arr(51, 21.5, 51, 15.8, "#2e7d32")
    box(32, 2, 38, 13, "#e8eef8", "#1f4e79")
    ax.text(51, 11.2, "Apical uptake  $k_{abs}=P_{eff}/(\delta\,ER)$", ha="center", fontsize=7.4)
    ax.text(51, 5.2, "$P_{eff}=P_0[1+E_{max}H]$", ha="center", fontsize=7.4)
    arr(70.8, 44, 79.2, 44, "#b00020"); ax.text(75, 46.3, "wash-out $k_{dil}$", ha="center", fontsize=6.6, color="#b00020")
    arr(70.8, 33, 79.2, 33, "#b00020"); ax.text(75, 35.3, "proteolysis", ha="center", fontsize=6.6, color="#b00020")
    box(79.5, 28, 19, 22, "#fbeaea", "#b00020"); ax.text(89, 39, "lost to lumen\nand protease", ha="center", va="center", fontsize=7.2)
    arr(70.8, 8.5, 79.2, 8.5); ax.text(75, 11.0, r"$\times S_{post}$", ha="center", fontsize=7)
    box(79.5, 2, 19, 13, "#f3e8fb", "#6a1b9a"); ax.text(89, 8.5, "Systemic $F$", ha="center", va="center", fontsize=7.8)
    arr(31.5, 8.5, 22.5, 8.5, "#555555")
    ax.text(1, 18, "Enhancer clearance $k_{dep}$ and\nepithelial recovery $t_{rev}$ set the\nduration of the exposure; the\ndamage-hours DH measure it", ha="left", va="center", fontsize=6.6, color="#444444")
    save(fig, "fig1_schematic")


def fig_arms():
    r = J("paired"); arms = ["R0_peptide_alone", "R1_reference_PE", "R2_reference_PE_efflux", "P1_dual_lock", "P2_no_efflux_module", "P3_no_colocalisation", "P4_no_protection", "P5_no_self_limiting_PE"]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})
    y = np.arange(len(arms))[::-1]
    for yi, a in zip(y, arms):
        s = r["arms"][a]; col = C["P"] if a.startswith("P1") else (C["R"] if a.startswith("R") else "#9ecae1")
        ax[0].plot([s["p05"] * 100, s["p95"] * 100], [yi, yi], color=col, lw=2.2, solid_capstyle="round"); ax[0].plot(s["median"] * 100, yi, "o", color="k", ms=4)
        ax[1].barh(yi, s["p_success"] * 100, xerr=1.96 * s["se"] * 100, color=col, capsize=2)
    ax[0].set_xscale("log"); ax[0].set_yticks(y); ax[0].set_yticklabels([LAB[a] for a in arms], fontsize=7)
    ax[0].axvline(4, ls="--", lw=0.9, color=C["Q"]); ax[0].text(4.2, len(arms) - 0.6, "4 %", color=C["Q"], fontsize=7)
    ax[0].set_xlabel("absolute bioavailability $F$ (%): median and 5-95 %"); ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    ax[1].set_yticks(y); ax[1].set_yticklabels([]); ax[1].set_xlabel("P($F\\geq4$ % and reversible), % (95 % MC error)"); ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig2_arms")


def fig_ablation():
    r = J("paired"); fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4))
    rem = ["P1_minus_delta", "P1_minus_k_dil", "P1_minus_mism", "P1_minus_k_rel", "P3_no_colocalisation", "P4_no_protection", "P5_no_self_limiting_PE", "P2_no_efflux_module", "P1_full_dose", "M0_matrix_only"]
    add = ["R1_plus_delta", "R1_plus_k_dil", "R1_plus_mism", "R1_plus_k_rel", "R1_plus_confinement_bundle", "R1_plus_prot", "R1_plus_k_dep", "R2_reference_PE_efflux", "R1_half_dose", "R1_plus_halfPE_efflux"]
    for a_, lst, key, ttl in ((ax[0], rem, "paired_vs_P1", "Removal from Dual-Lock P1"), (ax[1], add, "paired_vs_R1", "Addition to reference R1")):
        y = np.arange(len(lst))[::-1]
        for yi, k in zip(y, lst):
            d = r[key][k]["d_success"]; a_.errorbar(d[0] * 100, yi, xerr=[[ (d[0] - d[1]) * 100], [(d[2] - d[0]) * 100]], fmt="o", color=C["P"] if key == "paired_vs_P1" else C["G"], ms=4, capsize=2)
        a_.axvline(0, color="k", lw=0.7); a_.set_yticks(y); a_.set_yticklabels([LAB[k] for k in lst], fontsize=6.8); a_.set_title(ttl, fontsize=8)
        a_.set_xlabel("paired change in P(success), percentage points")
    fig.tight_layout(); save(fig, "fig3_ablation")


def fig_sobol():
    r = J("sobol"); fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.5), gridspec_kw={"width_ratios": [1.3, 1]})
    pp = r["P1_dual_lock"]["log10_F"]["params"]; keys = sorted(pp, key=lambda k: pp[k]["ST"], reverse=True)[:12][::-1]
    y = np.arange(len(keys)); h = 0.38
    ax[0].barh(y + h / 2, [pp[k]["ST"] for k in keys], h, xerr=[[pp[k]["ST"] - pp[k]["ST_ci"][0] for k in keys], [pp[k]["ST_ci"][1] - pp[k]["ST"] for k in keys]], color=C["P"], label="total-order $S_T$", capsize=1.5)
    ax[0].barh(y - h / 2, [pp[k]["S1"] for k in keys], h, xerr=[[pp[k]["S1"] - pp[k]["S1_ci"][0] for k in keys], [pp[k]["S1_ci"][1] - pp[k]["S1"] for k in keys]], color="#9ecae1", label="first-order $S_1$", capsize=1.5)
    nm = {"delta": "δ niche thickness", "P0": "$P_0$ baseline permeability", "EC50": "EC$_{50}$", "A_m": "$A_m$ contact area", "Emax": "$E_{max}$", "k_deg": "$k_{deg}$", "S_post": "$S_{post}$", "k_dil": "$k_{dil}$ wash-out",
          "h": "Hill $h$", "ER0": "$ER_0$", "sub": "substrate status", "IC50": "IC$_{50}$", "C_tox": "$C_{tox}$", "t_rev": "$t_{rev}$", "m": "release mismatch", "Pr": "protection", "k_dep": "$k_{dep}$", "k_rel": "$k_{rel}$"}
    ax[0].set_yticks(y); ax[0].set_yticklabels([nm.get(k, k) for k in keys], fontsize=7); ax[0].set_xlabel("Sobol index for log$_{10}F$ (P1)"); ax[0].legend(fontsize=6.5, frameon=False, loc="lower right")
    ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    gr = list(r["P1_dual_lock"]["log10_F"]["groups"].keys()); x = np.arange(len(gr)); w = 0.2
    for i, (arm, o, col) in enumerate((("P1_dual_lock", "log10_F", C["P"]), ("R1_reference_PE", "log10_F", C["R"]), ("P1_dual_lock", "success", "#9ecae1"), ("R1_reference_PE", "success", "#cfcfcf"))):
        v = [r[arm][o]["groups"][g]["total"] for g in gr]
        ax[1].bar(x + (i - 1.5) * w, v, w, color=col, label="%s, %s" % (arm[:2], "log$_{10}F$" if o == "log10_F" else "success"))
    ax[1].set_xticks(x); ax[1].set_xticklabels([g.replace(" and ", " &\n").replace(" ", "\n", 1) if len(g) > 14 else g for g in gr], fontsize=6.3)
    ax[1].set_ylabel("total-order group index"); ax[1].legend(fontsize=6, frameon=False); ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig4_sobol")


def fig_pareto():
    r = J("pareto"); pa = J("paired"); fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
    cols = {"R1 reference (SNAC-like)": C["R"], "R2 reference + efflux module": "#aaaaaa", "P1 Dual-Lock": C["P"], "P3 Dual-Lock without co-localisation": "#9ecae1", "P5 Dual-Lock without self-limiting PE": C["O"]}
    for name, rows in r["families"].items():
        x = [1 - s["p_safe"] for s in rows]; y = [s["p_F_ge_4"] for s in rows]
        ax[0].plot(np.array(x) * 100, np.array(y) * 100, "-o", ms=3, color=cols[name], label=name.replace("colocalisation", "co-localization"))
        if name in ("R1 reference (SNAC-like)", "P1 Dual-Lock"):
            for s in rows:
                if s["dose"] in (0.5, 1.0): ax[0].annotate("%g" % s["dose"], (100 * (1 - s["p_safe"]), 100 * s["p_F_ge_4"]), fontsize=6, xytext=(3, -8), textcoords="offset points", color=cols[name])
    ax[0].set_xlabel("P(not reversible), % (risk)"); ax[0].set_ylabel("P($F\\geq4$ %), % (benefit)"); ax[0].legend(fontsize=5.8, frameon=False, loc="lower right")
    ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    for k, c, lab in (("R1_reference_PE", C["R"], "R1"), ("P1_dual_lock", C["P"], "P1"), ("P5_no_self_limiting_PE", C["O"], "P5")):
        s = pa["samples"][k]; ax[1].scatter(np.maximum(s["DH"], 1e-4), np.array(s["F"]) * 100, s=3, alpha=0.25, color=c, label=lab, rasterized=True)
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].axvline(0.5, ls="--", color="k", lw=0.8); ax[1].axhline(4, ls="--", color=C["Q"], lw=0.8)
    ax[1].set_xlabel("damage-hours DH (h)"); ax[1].set_ylabel("$F$ (%)"); ax[1].set_ylim(1e-2, 100); ax[1].legend(fontsize=7, frameon=False, markerscale=3)
    ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig5_pareto")


def fig_fals():
    r = J("fals"); fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
    G = np.array(r["grid_delta_dose"]["p_success"]) * 100; B = np.array(r["region_B"], dtype=float)
    im = ax[0].imshow(G, origin="lower", aspect="auto", cmap="viridis", vmin=0, vmax=70)
    ax[0].contour(B, levels=[0.5], colors="w", linewidths=1.6)
    ax[0].set_xticks(range(len(r["deltas"]))); ax[0].set_xticklabels(r["deltas"], fontsize=6.5, rotation=45); ax[0].set_yticks(range(len(r["doses"]))); ax[0].set_yticklabels(r["doses"], fontsize=6.5)
    ax[0].set_xlabel("niche thickness δ (cm)"); ax[0].set_ylabel("PE dose (mmol)"); ax[0].set_title("P(success), %; white line: benefit region B", fontsize=7.5)
    for (i, j), v in np.ndenumerate(G):
        ax[0].text(j, i, "%.0f" % v, ha="center", va="center", fontsize=5.2, color="w" if v < 40 else "k")
    ax[0].text(-0.02, 1.08, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    H = np.array(r["grid_delta_kdil"]["p_success"]) * 100
    im2 = ax[1].imshow(H, origin="lower", aspect="auto", cmap="viridis", vmin=0, vmax=70)
    ax[1].set_xticks(range(len(r["deltas"]))); ax[1].set_xticklabels(r["deltas"], fontsize=6.5, rotation=45); ax[1].set_yticks(range(len(r["k_dils"]))); ax[1].set_yticklabels(r["k_dils"], fontsize=6.5)
    ax[1].set_xlabel("niche thickness δ (cm)"); ax[1].set_ylabel("niche wash-out $k_{dil}$ (h$^{-1}$)"); ax[1].set_title("P(success), %, at 0.5 mmol PE", fontsize=7.5)
    for (i, j), v in np.ndenumerate(H):
        ax[1].text(j, i, "%.0f" % v, ha="center", va="center", fontsize=5.2, color="w" if v < 40 else "k")
    ax[1].text(-0.02, 1.08, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.colorbar(im2, ax=ax, fraction=0.025, pad=0.02); save(fig, "fig6_falsification")


def fig_struct():
    r = J("struct")["variants"]; keys = [k for k in r if "P0 not recalibrated" not in k]
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1.5, 1]})
    x = np.arange(len(keys)); w = 0.38
    for i, k in enumerate(keys):
        v = r[k]
        if v.get("calibratable", True) is False:
            ax[0].text(i, 5, "reference\ncannot reach\n1.2 %", ha="center", fontsize=5.8, color=C["Q"]); continue
        ax[0].bar(i - w / 2, v["R1_reference_PE"]["p_success"] * 100, w, color=C["R"], label="R1" if i == 0 else None)
        ax[0].bar(i + w / 2, v["P1_dual_lock"]["p_success"] * 100, w, color=C["P"], label="P1" if i == 0 else None)
        ax[0].text(i, max(v["P1_dual_lock"]["p_success"], v["R1_reference_PE"]["p_success"]) * 100 + 2, "$P_0$×%.1f" % v["p0_over_base"], ha="center", fontsize=5.5)
    ax[0].set_xticks(x); ax[0].set_xticklabels([k.replace("spatial, dispersion ", "spatial\n").replace("lumped + stagnant layer", "lumped +\nlayer").replace("lumped (base)", "lumped\n(base)") for k in keys], fontsize=5.8)
    ax[0].set_ylabel("P(success), %"); ax[0].legend(fontsize=7, frameon=False); ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    ks = [k for k in keys if r[k].get("calibratable", True) is not False]
    ax[1].plot(range(len(ks)), [r[k]["P1_dual_lock"]["median"] * 100 for k in ks], "o-", color=C["P"], label="P1 median $F$")
    ax[1].plot(range(len(ks)), [r[k]["R1_reference_PE"]["median"] * 100 for k in ks], "s-", color=C["R"], label="R1 median $F$"); ax[1].axhline(4, ls="--", color=C["Q"], lw=0.8)
    ax[1].set_yscale("log"); ax[1].set_xticks(range(len(ks))); ax[1].set_xticklabels([k.replace("spatial, dispersion ", "×").replace("lumped + stagnant layer", "layer").replace("lumped (base)", "base") for k in ks], fontsize=6, rotation=45)
    ax[1].set_ylabel("median $F$ (%)"); ax[1].legend(fontsize=6.5, frameon=False); ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig7_structural")


def fig_conv_corr():
    cv = J("converge"); co = J("corr"); fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.8))
    for arm, c in (("R1_reference_PE", C["R"]), ("P1_dual_lock", C["P"])):
        rr = cv["running_" + arm]; n = np.array([d["n"] for d in rr]); p = np.array([d["p_success"] for d in rr]) * 100; se = np.sqrt(p / 100 * (1 - p / 100) / n) * 196
        ax[0].plot(n, p, "-o", ms=3, color=c, label=arm[:2]); ax[0].fill_between(n, p - se, p + se, color=c, alpha=0.2)
    ax[0].set_xscale("log"); ax[0].set_xlabel("number of draws"); ax[0].set_ylabel("P(success), %"); ax[0].legend(fontsize=7, frameon=False); ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    for arm, c in (("R1_reference_PE", C["R"]), ("P1_dual_lock", C["P"])):
        rr = cv["running_" + arm]; ax[1].plot([d["n"] for d in rr], [d["median_F"] * 100 for d in rr], "-o", ms=3, color=c)
    ax[1].set_xscale("log"); ax[1].set_xlabel("number of draws"); ax[1].set_ylabel("median $F$ (%)"); ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    rhos = sorted(co["sweep_rho_EC50_Ctox"], key=float)
    for arm, c in (("R1_reference_PE", C["R"]), ("P1_dual_lock", C["P"]), ("P5_no_self_limiting_PE", C["O"])):
        ax[2].plot([float(r_) for r_ in rhos], [co["sweep_rho_EC50_Ctox"][r_][arm]["p_success"] * 100 for r_ in rhos], "-o", ms=3, color=c, label=arm[:2])
    ax[2].set_xlabel("ρ(EC$_{50}$, C$_{tox}$) (Gaussian copula)"); ax[2].set_ylabel("P(success), %"); ax[2].legend(fontsize=7, frameon=False); ax[2].text(-0.02, 1.04, "c", transform=ax[2].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig8_convergence_correlation")


def fig_plaus():
    pa = J("paired"); ab = J("abc"); va = J("valid"); fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
    bins = np.logspace(-3, 0, 60)
    for k, c, lab in (("R1_reference_PE", C["R"], "R1 reference"), ("P1_dual_lock", C["P"], "P1 Dual-Lock")):
        ax[0].hist(np.maximum(pa["samples"][k]["F"], 1e-3) * 100 / 100, bins=bins, histtype="stepfilled", alpha=0.45, color=c, label=lab, density=True)
    ax[0].axvspan(0.004, 0.010, color=C["G"], alpha=0.18, label="oral semaglutide-like (0.4-1 %)"); ax[0].axvline(0.008, color=C["G"], lw=1.0)
    ax[0].axvline(0.06, color=C["Q"], lw=1.0, ls=":", label="MEDI7219 in dogs (~6 %)"); ax[0].axvline(0.04, color="k", lw=0.8, ls="--", label="4 % target")
    ax[0].set_xscale("log"); ax[0].set_xlabel("absolute bioavailability $F$"); ax[0].set_ylabel("density (log-$F$ bins)"); ax[0].legend(fontsize=5.8, frameon=False, loc="upper left"); ax[0].text(-0.02, 1.04, "a", transform=ax[0].transAxes, weight="bold", fontsize=10)
    tags = [k for k in ab if k not in ("n",)]
    vals = [ab[k]["P1"]["p_success"] * 100 for k in tags]; ns = [ab[k]["n_accepted"] for k in tags]
    ax[1].barh(range(len(tags))[::-1], vals, color=C["P"])
    for yi, v, n_ in zip(range(len(tags))[::-1], vals, ns): ax[1].text(v + 1, yi, "%.0f %% (n=%s)" % (v, "{:,}".format(n_)), va="center", fontsize=6.3)
    ax[1].set_yticks(range(len(tags))[::-1]); ax[1].set_yticklabels([t.replace(" (oral-semaglutide-like)", "") for t in tags], fontsize=6.5); ax[1].set_xlim(0, 85)
    ax[1].set_xlabel("P1 P(success) among accepted draws, %"); ax[1].text(-0.02, 1.04, "b", transform=ax[1].transAxes, weight="bold", fontsize=10)
    fig.tight_layout(); save(fig, "fig9_plausibility")


def fig_si_numerics():
    cv = J("converge"); fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.6))
    dts = sorted(cv["dt"], key=float); ax[0].plot([float(d) for d in dts], [cv["dt"][d]["median_F"] * 100 for d in dts], "o-", color=C["P"]); ax[0].set_xscale("log"); ax[0].set_xlabel("time step (h)"); ax[0].set_ylabel("P1 median $F$ (%)")
    for arm, c in (("R1_reference_PE", C["R"]), ("P1_dual_lock", C["P"])):
        hs = sorted(cv["horizon"][arm], key=float); ax[1].plot([float(h) for h in hs], [cv["horizon"][arm][h]["median_F"] * 100 for h in hs], "o-", color=c, label=arm[:2])
    ax[1].set_xlabel("horizon (h)"); ax[1].set_ylabel("median $F$ (%)"); ax[1].legend(fontsize=7, frameon=False)
    fig.tight_layout(); save(fig, "figS1_numerics")


if __name__ == "__main__":
    only = sys.argv[2:]
    fns = {"schematic": fig_schematic, "arms": fig_arms, "ablation": fig_ablation, "sobol": fig_sobol, "pareto": fig_pareto, "fals": fig_fals, "struct": fig_struct,
           "conv": fig_conv_corr, "plaus": fig_plaus, "numerics": fig_si_numerics}
    for k, f in fns.items():
        if only and k not in only: continue
        try: f()
        except FileNotFoundError as e: print("skip", k, e)
