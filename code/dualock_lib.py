# -*- coding: utf-8 -*-
"""dualock_lib.py -- extended machinery for the Dual-Lock oral-peptide Monte Carlo study (builds on oral_peptide_mc.py, which is unchanged).

* inputs_from_u : maps uniform variates (independent, Sobol or Gaussian-copula) to the priors of oral_peptide_mc.PRIORS / ARMS, so that every arm
                  can be run on the SAME random numbers (common random numbers -> paired contrasts).
* run_u         : lumped model on a U matrix, optionally with a serial stagnant-layer (mucus / unstirred-layer) resistance.
* integrate_pde : alternative structural model, 1-D diffusion-reaction across the niche (implicit finite volumes), same priors.
* sobol         : Saltelli / Jansen estimators (first- and total-order, single parameters and groups) with bootstrap intervals.
"""
import numpy as np
from scipy.special import ndtri, ndtr
from scipy.stats import beta as beta_dist, qmc, rankdata

import oral_peptide_mc as M

DIMS = ["P0", "k_deg", "S_post", "ER0", "sub", "A_m", "EC50", "Emax", "h", "C_tox", "t_rev", "IC50",
        "delta", "k_dil", "m", "Pr", "k_dep", "k_rel"]
IDX = {k: i for i, k in enumerate(DIMS)}
GROUPS = {"peptide and epithelium": ["P0", "k_deg", "S_post", "ER0", "sub", "A_m"],
          "enhancer pharmacology": ["EC50", "Emax", "h", "C_tox", "t_rev"],
          "efflux modulator": ["IC50"],
          "formulation architecture": ["delta", "k_dil", "m", "Pr", "k_dep", "k_rel"]}
SEC_PER_H = 3600.0


def arm_spec(name, **over):
    a = dict(M.ARMS[name]); a.update(over); return a


# ------------------------------------------------------------------------------------------------ inputs
def _ln(u, med, sig):
    return med * np.exp(sig * ndtri(np.clip(u, 1e-12, 1 - 1e-12)))


def _beta(u, mean, sd):
    k = mean * (1 - mean) / sd ** 2 - 1.0
    return beta_dist.ppf(np.clip(u, 1e-12, 1 - 1e-12), mean * k, (1 - mean) * k)


def inputs_from_u(U, a, p0, p_sub=None, emax_mult=1.0, sig_scale=1.0, extra=None):
    """U: (n, 18) uniforms in the column order DIMS. a: arm dict (medians of formulation parameters). Returns the input dict used by M.integrate."""
    s = sig_scale; P = M.PRIORS; g = {}
    c = lambda k: U[:, IDX[k]]
    g["P0"] = _ln(c("P0"), p0, P["P0"]["sig"] * s)
    g["k_deg"] = _ln(c("k_deg"), P["k_deg"]["med"], P["k_deg"]["sig"] * s)
    g["S_post"] = _beta(c("S_post"), P["S_post"]["mean"], P["S_post"]["sd"] * s)
    ps = P["p_sub"]["val"] if p_sub is None else p_sub
    sub = c("sub") < ps
    er0 = _ln(c("ER0"), P["ER0"]["med"], P["ER0"]["sig"] * s)
    g["ER0"] = np.where(sub, np.maximum(er0, 1.0), 1.0)
    g["A_m"] = _ln(c("A_m"), P["A_m"]["med"], P["A_m"]["sig"] * s)
    g["EC50"] = _ln(c("EC50"), P["EC50"]["med"], P["EC50"]["sig"] * s)
    g["Emax"] = _ln(c("Emax"), P["Emax"]["med"] * emax_mult, P["Emax"]["sig"] * s)
    g["h"] = _ln(c("h"), P["h"]["med"], P["h"]["sig"] * s)
    g["C_tox"] = _ln(c("C_tox"), P["C_tox"]["med"], P["C_tox"]["sig"] * s)
    g["t_rev"] = _ln(c("t_rev"), P["t_rev"]["med"], P["t_rev"]["sig"] * s)
    g["IC50"] = _ln(c("IC50"), P["IC50"]["med"], P["IC50"]["sig"] * s)
    g["delta"] = _ln(c("delta"), a["delta"], 0.30 * s)
    g["k_dil"] = _ln(c("k_dil"), a["k_dil"], 0.30 * s)
    g["m"] = _ln(c("m"), 1.0, a["mism"] * s)
    g["Pr"] = _beta(c("Pr"), max(a["prot"], 1e-3), 0.10 * s) if a["prot"] > 0 else np.zeros(len(U))
    g["k_dep"] = _ln(c("k_dep"), a["k_dep"], 0.30 * s)
    g["k_rel"] = _ln(c("k_rel"), a["k_rel"], 0.20 * s)
    if extra:
        g.update(extra)
    return g


def uniforms(n, seed, corr=None):
    """Independent uniforms, or a Gaussian copula with the given correlation dict {(name_i, name_j): rho} on the normal scores."""
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((n, len(DIMS)))
    if corr:
        R = np.eye(len(DIMS))
        for (i, j), r in corr.items():
            R[IDX[i], IDX[j]] = R[IDX[j], IDX[i]] = r
        w, V = np.linalg.eigh(R)
        if w.min() < 1e-8:                                   # not positive definite: nearest correlation matrix by eigenvalue clipping
            R = (V * np.clip(w, 1e-6, None)) @ V.T; d = np.sqrt(np.diag(R)); R = R / np.outer(d, d)
        Z = Z @ np.linalg.cholesky(R).T
    return ndtr(Z)


def sobol_uniforms(n_pow2, seed):
    sam = qmc.Sobol(d=2 * len(DIMS), scramble=True, seed=seed).random_base2(n_pow2)
    return sam[:, :len(DIMS)], sam[:, len(DIMS):]


# ------------------------------------------------------------------------------------------------ lumped model (+ optional stagnant layer)
SOL = dict(D_mg=100.0, D_sig=0.35, S_mg_ml=20.0, S_sig=0.8)   # peptide dose (50-200 mg) and solubility (grade D)


def run_u(U, spec, p0, pe=None, inh=None, p_sub=None, emax_mult=1.0, sig_scale=1.0, layer=None, dt=M.DT, t_end=M.T_END, track=False, sol=None, mods=None):
    """Lumped model on uniforms U. layer=(U_layer (n,2)) adds a serial stagnant layer: P_eff -> 1/(1/P_eff + h_s/D_m) with h_s, D_m lognormal.
    sol=(U_sol (n,2)) adds a peptide-solubility limit on the dissolved amount in the niche (P <= S V / D). mods: multipliers applied to named inputs (e.g. {'C_tox': 0.3})."""
    extra = None
    if layer is not None:
        hs = _ln(layer[:, 0], LAYER["h_s_cm"], LAYER["h_s_sig"])
        Dm = _ln(layer[:, 1], LAYER["D_cm2_s"], LAYER["D_sig"]) * SEC_PER_H
        extra = {"Ps": Dm / hs}
    if sol is not None:
        extra = dict(extra or {}); extra["D_mg"] = _ln(sol[:, 0], SOL["D_mg"], SOL["D_sig"]); extra["S_sol"] = _ln(sol[:, 1], SOL["S_mg_ml"], SOL["S_sig"])
    g = inputs_from_u(U, spec, p0, p_sub, emax_mult, sig_scale, extra)
    for k, v in (mods or {}).items():
        g[k] = g[k] * v
    return g, integrate_serial(g, spec, pe, inh, dt, t_end, track)


LAYER = dict(h_s_cm=0.015, h_s_sig=0.7, D_cm2_s=1.0e-6, D_sig=0.5)   # stagnant (mucus / unstirred) layer: Atuma 2001 (20-830 um); Cu & Saltzman 2009


def integrate_serial(g, spec, pe=None, inh=None, dt=M.DT, t_end=M.T_END, track=False):
    """Same equations as oral_peptide_mc.integrate, except that when g contains 'Ps' the permeability is P_eff,s = 1/(1/P_eff + 1/Ps)."""
    if "Ps" not in g and "S_sol" not in g:
        return M.integrate(g, "P1_dual_lock", dt=dt, t_end=t_end, pe=spec["pe"] if pe is None else pe, inh=spec["inh"] if inh is None else inh, track=track)
    n = len(g["P0"]); V = g["A_m"] * g["delta"]
    pe_dose = spec["pe"] if pe is None else pe; inh_dose = spec["inh"] if inh is None else inh
    Pt = np.ones(n); P = np.zeros(n); Et = np.full(n, pe_dose); E = np.zeros(n); It = np.full(n, inh_dose); I = np.zeros(n)
    Abs = np.zeros(n); Deg = np.zeros(n); Dil = np.zeros(n); DH = np.zeros(n); t_last = np.zeros(n)
    kdeg_eff = g["k_deg"] * (1.0 - g["Pr"])
    Pmax = g["S_sol"] * V / g["D_mg"] if "S_sol" in g else None
    for i in range(int(round(t_end / dt))):
        t = i * dt
        Ce = 1000.0 * E / V; Ci = 1000.0 * I / V
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(Ce > 0, (Ce / g["EC50"]) ** g["h"], 0.0); H = x / (1.0 + x)
            y = np.where(Ce > 0, (Ce / g["C_tox"]) ** 3.0, 0.0); Hd = y / (1.0 + y)
        Peff = g["P0"] * (1.0 + g["Emax"] * H)
        if "Ps" in g:
            Peff = 1.0 / (1.0 / Peff + 1.0 / g["Ps"])
        ER = 1.0 + (g["ER0"] - 1.0) / (1.0 + Ci / g["IC50"])
        kabs = Peff / (g["delta"] * ER)
        rP = g["k_rel"] * Pt
        if "S_sol" in g:
            rP = rP * np.clip(1.0 - P / Pmax, 0.0, 1.0)
        rE = g["k_rel"] * g["m"] * Et; rI = g["k_rel"] * g["m"] * It
        a_flux = kabs * P; d_flux = kdeg_eff * P; w_flux = g["k_dil"] * P
        Pt = Pt - rP * dt; P = P + (rP - a_flux - d_flux - w_flux) * dt
        Abs += a_flux * dt; Deg += d_flux * dt; Dil += w_flux * dt
        Et = Et - rE * dt; E = E + (rE - (g["k_dil"] + g["k_dep"]) * E) * dt
        It = It - rI * dt; I = I + (rI - g["k_dil"] * I) * dt
        DH += Hd * dt; t_last = np.where(Ce > 0.1 * g["EC50"], t + dt, t_last)
    F = Abs * g["S_post"]
    safe = (DH <= M.DH_MAX) & ((t_last + g["t_rev"]) <= M.T_REVERSE_MAX)
    return dict(F=F, safe=safe, DH=DH, t_last=t_last, Abs=Abs, Deg=Deg, Dil=Dil, rem=Pt + P, closure=Pt + P + Abs + Deg + Dil)


def success(o):
    return (o["F"] >= M.F_TARGET) & o["safe"]


def summ(o):
    F = o["F"]; ok = success(o)
    return dict(median=float(np.median(F)), p05=float(np.percentile(F, 5)), p95=float(np.percentile(F, 95)), p_F_ge_4=float(np.mean(F >= M.F_TARGET)),
                p_safe=float(np.mean(o["safe"])), p_success=float(ok.mean()), se=float(np.sqrt(ok.mean() * (1 - ok.mean()) / len(F))),
                median_DH=float(np.median(o["DH"])), p90_DH=float(np.percentile(o["DH"], 90)), median_trev=float(np.median(o["t_last"])))


# ------------------------------------------------------------------------------------------------ spatial (diffusion-reaction) model
DIFF = dict(P_cm2_s=1.0e-6, P_sig=0.5, E_cm2_s=4.0e-6, E_sig=0.5, I_cm2_s=2.0e-6, I_sig=0.5)


def _thomas(lo, di, up, rhs):
    """Vectorised tridiagonal solve over axis 0 (draws); arrays (n, nx); lo[:,0] and up[:,-1] unused."""
    n, nx = di.shape
    cp = np.empty_like(di); dp = np.empty_like(di)
    cp[:, 0] = up[:, 0] / di[:, 0]; dp[:, 0] = rhs[:, 0] / di[:, 0]
    for i in range(1, nx):
        den = di[:, i] - lo[:, i] * cp[:, i - 1]
        cp[:, i] = up[:, i] / den if i < nx - 1 else 0.0
        dp[:, i] = (rhs[:, i] - lo[:, i] * dp[:, i - 1]) / den
    x = np.empty_like(di); x[:, -1] = dp[:, -1]
    for i in range(nx - 2, -1, -1):
        x[:, i] = dp[:, i] - cp[:, i] * x[:, i + 1]
    return x


def integrate_pde(g, spec, pe=None, inh=None, nx=24, dt=0.02, t_end=M.T_END, mix=1.0, Dfac=None, return_profile=False):
    """1-D diffusion-reaction across the niche (x = 0 matrix / lumen side, x = delta epithelial surface).
    Peptide: diffusion D_P*mix; first-order proteolysis k_deg(1-Pr); wash-out at x=0 through a mass-transfer coefficient h0 = k_dil*delta (cm/h)
    (so that the well-mixed limit is the lumped model); release into the first cell; absorption at x=delta with permeability P_eff/ER, where P_eff and ER
    are evaluated from the enhancer / modulator concentration in the epithelial cell. Enhancer and modulator: same geometry, D_E, D_I; clearance k_dep
    (enhancer only). Backward Euler (unconditionally stable, exactly mass-conservative); enhancer/modulator advanced before the peptide at each step."""
    n = len(g["P0"]); pe_dose = spec["pe"] if pe is None else pe; inh_dose = spec["inh"] if inh is None else inh
    dxx = g["delta"] / nx                                  # (n,)
    DP = g["D_P"] * mix; DE = g["D_E"] * mix; DI = g["D_I"] * mix
    h0 = g["k_dil"] * g["delta"]
    kd = g["k_deg"] * (1.0 - g["Pr"])
    cP = np.zeros((n, nx)); cE = np.zeros((n, nx)); cI = np.zeros((n, nx))
    Pt = np.ones(n); Et = np.full(n, pe_dose); It = np.full(n, inh_dose)
    Abs = np.zeros(n); Deg = np.zeros(n); Dil = np.zeros(n); DH = np.zeros(n); t_last = np.zeros(n)
    A = g["A_m"]

    def system(D, loss0, lossN, react):
        a = D / dxx ** 2
        lo = np.repeat((-dt * a)[:, None], nx, 1); up = lo.copy()
        di = np.repeat((1 + dt * 2 * a)[:, None], nx, 1)
        di[:, 0] = 1 + dt * (a + loss0)
        di[:, -1] = 1 + dt * (a + lossN)
        di += dt * react[:, None]
        return lo, di, up

    heffE = 1.0 / (1.0 / h0 + dxx / (2 * DE)); heffI = 1.0 / (1.0 / h0 + dxx / (2 * DI)); heffP = 1.0 / (1.0 / h0 + dxx / (2 * DP))
    loE, diE, upE = system(DE, heffE / dxx, 0.0, g["k_dep"])
    loI, diI, upI = system(DI, heffI / dxx, 0.0, np.zeros(n))
    steps = int(round(t_end / dt))
    for i in range(steps):
        t = i * dt
        rE = g["k_rel"] * g["m"] * Et; rI = g["k_rel"] * g["m"] * It; rP = g["k_rel"] * Pt
        sE = np.zeros((n, nx)); sE[:, 0] = rE / (A * dxx); sI = np.zeros((n, nx)); sI[:, 0] = rI / (A * dxx); sP = np.zeros((n, nx)); sP[:, 0] = rP / (A * dxx)
        Et = Et - rE * dt; It = It - rI * dt; Pt = Pt - rP * dt
        cE = _thomas(loE, diE, upE, cE + dt * sE); cI = _thomas(loI, diI, upI, cI + dt * sI)
        Ce = 1000.0 * cE[:, -1]; Ci = 1000.0 * cI[:, -1]
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(Ce > 0, (Ce / g["EC50"]) ** g["h"], 0.0); Hh = x / (1.0 + x)
            y = np.where(Ce > 0, (Ce / g["C_tox"]) ** 3.0, 0.0); Hd = y / (1.0 + y)
        Peff = g["P0"] * (1.0 + g["Emax"] * Hh)
        if "Ps" in g:
            Peff = 1.0 / (1.0 / Peff + 1.0 / g["Ps"])
        ER = 1.0 + (g["ER0"] - 1.0) / (1.0 + Ci / g["IC50"])
        kp = Peff / ER
        keff = 1.0 / (1.0 / kp + dxx / (2 * DP))
        loP, diP, upP = system(DP, heffP / dxx, keff / dxx, kd)
        cP = _thomas(loP, diP, upP, cP + dt * sP)
        Abs += dt * A * keff * cP[:, -1]
        Deg += dt * A * (kd * cP.sum(1) * dxx)
        Dil += dt * A * heffP * cP[:, 0]
        DH += Hd * dt; t_last = np.where(Ce > 0.1 * g["EC50"], t + dt, t_last)
    F = Abs * g["S_post"]
    safe = (DH <= M.DH_MAX) & ((t_last + g["t_rev"]) <= M.T_REVERSE_MAX)
    inv = A * dxx * cP.sum(1)
    out = dict(F=F, safe=safe, DH=DH, t_last=t_last, Abs=Abs, Deg=Deg, Dil=Dil, rem=Pt + inv, closure=Pt + inv + Abs + Deg + Dil)
    if return_profile:
        out["cP"] = cP; out["cE"] = cE
    return out


def run_pde(U, spec, p0, UD, mix=1.0, pe=None, inh=None, p_sub=None, nx=24, dt=0.02, layer=None, **kw):
    """UD: (n,3) uniforms for D_P, D_E, D_I."""
    extra = {"D_P": _ln(UD[:, 0], DIFF["P_cm2_s"], DIFF["P_sig"]) * SEC_PER_H, "D_E": _ln(UD[:, 1], DIFF["E_cm2_s"], DIFF["E_sig"]) * SEC_PER_H,
             "D_I": _ln(UD[:, 2], DIFF["I_cm2_s"], DIFF["I_sig"]) * SEC_PER_H}
    g = inputs_from_u(U, spec, p0, p_sub, 1.0, 1.0, extra)
    return g, integrate_pde(g, spec, pe, inh, nx=nx, dt=dt, mix=mix, **kw)


def calibrate_generic(fun, target=M.F_BASELINE_MID, iters=22, lo=1e-6, hi=2.0):
    """Bisection on the baseline permeability median so that median F of fun(p0) equals the target (same anchor as the lumped model)."""
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        if np.median(fun(mid)["F"]) > target:
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


# ------------------------------------------------------------------------------------------------ Sobol
def sobol(fun, n_pow2, seed, groups=None, nboot=200, rng_seed=0):
    """Saltelli (2010) first-order and Jansen (1999) total-order indices for the vector-valued output returned by fun(U) -> {name: y}.
    Returns {output: {param: {S1, ST, S1_ci, ST_ci}}} and, for groups, closed (first-order group) and total-group indices."""
    A, B = sobol_uniforms(n_pow2, seed)
    n = len(A)
    yA = fun(A); yB = fun(B)
    names = list(yA.keys())
    AB = {}
    for i, k in enumerate(DIMS):
        C = A.copy(); C[:, i] = B[:, i]; AB[k] = fun(C)
    GAB = {}
    if groups:
        for gname, members in groups.items():
            C = A.copy()
            for k in members:
                C[:, IDX[k]] = B[:, IDX[k]]
            GAB[gname] = fun(C)
    rng = np.random.default_rng(rng_seed)
    boots = [rng.integers(0, n, n) for _ in range(nboot)]

    def est(a, b, ab, sel):
        a, b, ab = a[sel], b[sel], ab[sel]
        var = np.var(np.concatenate([a, b]), ddof=1)
        return float(np.mean(b * (ab - a)) / var), float(0.5 * np.mean((a - ab) ** 2) / var)

    res = {}
    for o in names:
        res[o] = {"params": {}, "groups": {}}
        for k in DIMS:
            s1, st = est(yA[o], yB[o], AB[k][o], slice(None))
            bs = np.array([est(yA[o], yB[o], AB[k][o], b) for b in boots])
            res[o]["params"][k] = dict(S1=s1, ST=st, S1_ci=[float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))],
                                       ST_ci=[float(np.percentile(bs[:, 1], 2.5)), float(np.percentile(bs[:, 1], 97.5))])
        for gname in GAB:
            s1, st = est(yA[o], yB[o], GAB[gname][o], slice(None))
            bs = np.array([est(yA[o], yB[o], GAB[gname][o], b) for b in boots])
            res[o]["groups"][gname] = dict(closed=s1, total=st, closed_ci=[float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))],
                                           total_ci=[float(np.percentile(bs[:, 1], 2.5)), float(np.percentile(bs[:, 1], 97.5))])
    return res
