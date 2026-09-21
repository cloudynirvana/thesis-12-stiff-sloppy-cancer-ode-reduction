#!/usr/bin/env python3
"""Stiff-sloppy spectra and one MBAM-style reduction for a 10-parameter cancer-state ODE.

Seed 20260921. The Fisher matrices are noise-free (Gaussian information of the mean).
One noisy draw is used only for the observed-Hessian check and the flat-canyon scan.
Research only. Not a medical device, not a dose, not a cell-line fit.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260921
RNG = np.random.default_rng(SEED)

NAMES = ["r", "a", "aq", "K", "q", "f", "c", "gp", "gs", "h"]
# Proliferation, apoptosis of P, apoptosis of Q, carrying scale,
# basal P→Q, Q→P return, apoptotic clearance, stress production,
# stress clearance, stress gain on the P→Q switch.
THETA0 = np.array([0.42, 0.11, 0.025, 5.5, 0.22, 0.06, 0.85, 0.40, 2.40, 1.10])
Y0 = np.array([0.35, 0.04, 0.0, 0.02])  # P, Q, A, S
T_END = 40.0
T_OBS = np.array([2.0, 4.0, 7.0, 11.0, 16.0, 22.0, 30.0, 40.0])
T_CLAIM_EARLY = 11.0
T_CLAIM_LATE = 40.0

# Declared before the spectra were read as a scientific conclusion.
# Numerical zero: 1e-8 times the leading eigenvalue.
# Sloppy cut: 1e-3 times the leading eigenvalue (practical sloppiness).
NUM_CUT = 1e-8
SLOPPY_CUT = 1e-3
FD_STEP = 1e-4
FLOW_SPEED = 0.22
FLOW_MAX_STEPS = 36
FLOW_FACTOR_STOP = 30.0
FLOW_RIDGE = 0.02

SIGMA = {"B": 0.05, "L": 0.08}
N_REP = {"B": 4, "L": 2}
A_FLOOR = 1e-4

# Invented rebound coefficient in the bad vector field. Not estimated.
BAD_RHO = 0.35


def unpack(th):
    r, a, aq, K, q, f, c, gp, gs, h = [float(v) for v in th]
    return r, a, aq, K, q, f, c, gp, gs, h


def rhs_full(t, y, th):
    P, Q, A, S = y
    r, a, aq, K, q, f, c, gp, gs, h = unpack(th)
    B = P + Q
    switch = q * (1.0 + h * S)
    dP = r * P * (1.0 - B / K) - switch * P + f * Q - a * P
    dQ = switch * P - f * Q - aq * Q
    dA = a * P + aq * Q - c * A
    dS = gp * B - gs * S
    return [dP, dQ, dA, dS]


def rhs_qss(t, y, th):
    """Quasi-steady stress. kappa = h * gp / gs is one parameter.

    th order: r, a, aq, K, q, f, c, kappa
    States: P, Q, A. S is slaved, S = (gp/gs) * B, and h*S = kappa * B.
    Clearance c remains, so map L can still be formed. Map B does not use A or c.
    """
    P, Q, A = y
    r, a, aq, K, q, f, c, kappa = [float(v) for v in th]
    B = P + Q
    switch = q * (1.0 + kappa * B)
    dP = r * P * (1.0 - B / K) - switch * P + f * Q - a * P
    dQ = switch * P - f * Q - aq * Q
    dA = a * P + aq * Q - c * A
    return [dP, dQ, dA]


def rhs_bad(t, y, th):
    """Bad reduction: apoptotic clearance is rewritten as a rebound drive on P.

    The original map-B Fisher column for c is empty. This vector field spends that
    free direction on a new term rho * c * A inside dP/dt and then treats the term
    as a discovered mechanism.
    """
    P, Q, A, S = y
    r, a, aq, K, q, f, c, gp, gs, h = unpack(th)
    B = P + Q
    switch = q * (1.0 + h * S)
    dP = r * P * (1.0 - B / K) - switch * P + f * Q - a * P + BAD_RHO * c * A
    dQ = switch * P - f * Q - aq * Q
    dA = a * P + aq * Q - c * A
    dS = gp * B - gs * S
    return [dP, dQ, dA, dS]


def _integrate(fun, y0, th, t_eval):
    sol = solve_ivp(
        fun,
        (0.0, float(t_eval[-1])),
        y0,
        t_eval=np.asarray(t_eval, dtype=float),
        args=(th,),
        method="LSODA",
        rtol=1e-8,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


def trajectory(th, t_eval=None, kind="full"):
    if t_eval is None:
        t_eval = np.linspace(0.0, T_END, 401)
    th = np.asarray(th, dtype=float)
    if kind == "full":
        return _integrate(rhs_full, Y0, th, t_eval)
    if kind == "bad":
        return _integrate(rhs_bad, Y0, th, t_eval)
    if kind == "qss":
        yq = _integrate(rhs_qss, Y0[:3], th, t_eval)
        return yq
    raise ValueError(kind)


def observe_from_state(state, kind, amap):
    """state columns are times. kind 'full'/'bad' use P,Q,A; 'qss' uses P,Q,A without S."""
    P = state[0]
    Q = state[1]
    A = state[2]
    if amap == "B":
        return np.log(np.maximum(P + Q, 1e-12))
    if amap == "L":
        return np.concatenate(
            [
                np.log(np.maximum(P, 1e-12)),
                np.log(np.maximum(Q, 1e-12)),
                np.log(np.maximum(A, 0.0) + A_FLOOR),
            ]
        )
    raise ValueError(amap)


def predict(th, amap, kind="full"):
    state = trajectory(np.asarray(th, dtype=float), T_OBS, kind=kind)
    return observe_from_state(state, kind, amap)


def burden_at(th, t, kind="full"):
    state = trajectory(np.asarray(th, dtype=float), [t], kind=kind)
    return float(state[0, 0] + state[1, 0])


def claim_log_ratio(th, kind="full"):
    b1 = burden_at(th, T_CLAIM_EARLY, kind=kind)
    b2 = burden_at(th, T_CLAIM_LATE, kind=kind)
    return float(np.log(b2) - np.log(b1))


def kappa_of(th):
    r, a, aq, K, q, f, c, gp, gs, h = unpack(th)
    return h * gp / gs


def qss_theta(th):
    r, a, aq, K, q, f, c, gp, gs, h = unpack(th)
    return np.array([r, a, aq, K, q, f, c, h * gp / gs])


def fisher(th, amap):
    """Gaussian Fisher matrix in log-parameters. Noise-free sensitivities."""
    th = np.asarray(th, dtype=float)
    mu = predict(th, amap, kind="full")
    n = len(th)
    sens = np.zeros((mu.size, n))
    for j in range(n):
        step = FD_STEP
        up = th.copy()
        up[j] = th[j] * np.exp(step)
        dn = th.copy()
        dn[j] = th[j] * np.exp(-step)
        # derivative w.r.t. log theta_j
        sens[:, j] = (predict(up, amap) - predict(dn, amap)) / (2.0 * step)
    weight = N_REP[amap] / (SIGMA[amap] ** 2)
    fim = weight * (sens.T @ sens)
    return 0.5 * (fim + fim.T), sens, mu


def spectrum_report(fim):
    evals, evecs = np.linalg.eigh(fim)
    # eigh is ascending
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]
    lam_max = float(evals[0])
    ratios = evals / lam_max
    numerical_rank = int(np.sum(ratios > NUM_CUT))
    practical_rank = int(np.sum(ratios > SLOPPY_CUT))
    sloppy_idx = np.where(ratios <= SLOPPY_CUT)[0]
    retained = evals[ratios > SLOPPY_CUT]
    if retained.size >= 2 and retained[-1] > 0:
        cond = float(retained[0] / retained[-1])
    elif retained.size == 1:
        cond = 1.0
    else:
        cond = None
    participation = (evecs[:, sloppy_idx] ** 2).sum(axis=1) if sloppy_idx.size else np.zeros(len(NAMES))
    return {
        "eigenvalues_desc": evals.tolist(),
        "eigenvectors_desc_columns": evecs.tolist(),
        "lambda_max": lam_max,
        "ratios_to_max": ratios.tolist(),
        "numerical_rank": numerical_rank,
        "practical_rank": practical_rank,
        "n_param": len(NAMES),
        "condition_of_stiff_block": cond,
        "sloppy_participation": participation.tolist(),
        "sloppy_count": int(sloppy_idx.size),
    }


def claim_gate(th, amap, spec):
    """Fraction of ||d g / d log theta||^2 that lies in the sloppy subspace."""
    th = np.asarray(th, dtype=float)
    g0 = claim_log_ratio(th, kind="full")
    grad = np.zeros(len(th))
    for j in range(len(th)):
        up = th.copy()
        up[j] = th[j] * np.exp(FD_STEP)
        dn = th.copy()
        dn[j] = th[j] * np.exp(-FD_STEP)
        grad[j] = (claim_log_ratio(up) - claim_log_ratio(dn)) / (2.0 * FD_STEP)
    evecs = np.array(spec["eigenvectors_desc_columns"])
    ratios = np.array(spec["ratios_to_max"])
    stiff = evecs[:, ratios > SLOPPY_CUT]
    sloppy = evecs[:, ratios <= SLOPPY_CUT]
    gnorm2 = float(grad @ grad)
    stiff_energy = float(np.sum((stiff.T @ grad) ** 2)) if stiff.size else 0.0
    sloppy_energy = float(np.sum((sloppy.T @ grad) ** 2)) if sloppy.size else 0.0
    return {
        "g_full": g0,
        "log_grad": grad.tolist(),
        "grad_norm2": gnorm2,
        "fraction_in_sloppy_subspace": (sloppy_energy / gnorm2) if gnorm2 > 0 else None,
        "fraction_in_stiff_subspace": (stiff_energy / gnorm2) if gnorm2 > 0 else None,
        "gate_cut": SLOPPY_CUT,
        "gated": bool(gnorm2 > 0 and sloppy_energy / gnorm2 < 0.05),
    }


def cost_against_reference(th, amap, mu_ref, kind="full"):
    resid = predict(th, amap, kind=kind) - mu_ref
    weight = N_REP[amap] / (SIGMA[amap] ** 2)
    return float(weight * np.dot(resid, resid))


def eigenvector_flow(th0, amap, mu_ref):
    """Sloppiest-eigenvector flow with a least-squares correction off the ray.

    This is not the Christoffel geodesic. It is the approximation stated in the methods:
    step along the current smallest eigenvector, then refit the orthogonal complement
    to the noise-free observations of the full model.
    """
    phi = np.log(np.asarray(th0, dtype=float))
    prev = None
    rows = []
    for step in range(FLOW_MAX_STEPS + 1):
        th = np.exp(phi)
        fim, _, _ = fisher(th, amap)
        evals, evecs = np.linalg.eigh(fim)
        v = evecs[:, 0].copy()
        if prev is not None and float(np.dot(v, prev)) < 0.0:
            v = -v
        ratio = np.exp(phi - np.log(th0))
        chi = cost_against_reference(th, amap, mu_ref, kind="full")
        rows.append(
            {
                "step": step,
                "chi2": chi,
                "log10_smallest_eig": float(np.log10(max(float(evals[0]), 1e-30))),
                "smallest_eig": float(max(float(evals[0]), 0.0)),
                "parameter_ratio_to_start": ratio.tolist(),
                "eigenvector": v.tolist(),
                "kappa": kappa_of(th),
                "max_abs_log_ratio": float(np.max(np.abs(np.log(ratio)))),
            }
        )
        if step == FLOW_MAX_STEPS:
            break
        if np.max(ratio) > FLOW_FACTOR_STOP or np.min(ratio) < 1.0 / FLOW_FACTOR_STOP:
            break
        prev = v.copy()
        proposal = phi + FLOW_SPEED * v
        # Refit in the orthogonal complement so the step stays near the model manifold.
        basis = evecs[:, 1:]  # Euclidean complement of v; columns

        def residual(coef):
            trial = np.exp(proposal + basis @ coef)
            data = predict(trial, amap) - mu_ref
            # Ridge on the correction. Map B has fewer residuals than free coefficients.
            return np.concatenate([data, FLOW_RIDGE * coef])

        fit = least_squares(
            residual,
            np.zeros(basis.shape[1]),
            method="trf",
            xtol=1e-12,
            ftol=1e-12,
        )
        phi = proposal + basis @ fit.x
    return rows


def scan_parameter(th0, amap, mu_ref, name, factors):
    j = NAMES.index(name)
    out = []
    for fac in factors:
        th = np.array(th0, dtype=float)
        th[j] = th0[j] * fac
        out.append(
            {
                "factor": float(fac),
                "value": float(th[j]),
                "chi2": cost_against_reference(th, amap, mu_ref),
                "g": claim_log_ratio(th),
            }
        )
    return out


def observed_hessian(th0, amap, mu_noisy):
    """Finite-difference Hessian of chi-square/2? We store Hessian of the weighted RSS.

    C(phi) = weight * ||mu(phi) - mu_noisy||^2
    At a noise-free match, C's Hessian is 2 * Fisher if Fisher was built as weight * J^T J.
    We compare eigenvalues of 0.5 * H to the Fisher eigenvalues.
    """
    phi0 = np.log(np.asarray(th0, dtype=float))
    h = 2e-4
    n = len(phi0)

    def cost(phi):
        resid = predict(np.exp(phi), amap) - mu_noisy
        weight = N_REP[amap] / (SIGMA[amap] ** 2)
        return float(weight * np.dot(resid, resid))

    H = np.zeros((n, n))
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h
        for j in range(i, n):
            ej = np.zeros(n)
            ej[j] = h
            if i == j:
                val = cost(phi0 + ei) - 2.0 * cost(phi0) + cost(phi0 - ei)
                H[i, j] = val / (h ** 2)
            else:
                val = (
                    cost(phi0 + ei + ej)
                    - cost(phi0 + ei - ej)
                    - cost(phi0 - ei + ej)
                    + cost(phi0 - ei - ej)
                )
                H[i, j] = val / (4.0 * h * h)
                H[j, i] = H[i, j]
    H = 0.5 * (H + H.T)
    # Fisher in this file is weight * J^T J, which equals (1/2) Hessian of C when residual is 0.
    half = 0.5 * H
    evals = np.linalg.eigvalsh(half)
    evals = np.sort(evals)[::-1]
    return half, evals


def rel_err(a, b):
    return float(abs(a - b) / max(abs(b), 1e-12))


def main():
    th0 = THETA0.copy()
    mu = {}
    specs = {}
    fishers = {}
    gates = {}
    for amap in ("B", "L"):
        fim, sens, mu_ref = fisher(th0, amap)
        fishers[amap] = fim
        mu[amap] = mu_ref
        spec = spectrum_report(fim)
        specs[amap] = spec
        gates[amap] = claim_gate(th0, amap, spec)
        print(
            f"map {amap}: rank {spec['numerical_rank']}/{spec['n_param']} "
            f"practical {spec['practical_rank']} cond {spec['condition_of_stiff_block']}"
        )
        print("  eigenvalues", np.array2string(np.array(spec["eigenvalues_desc"]), precision=4))
        evecs = np.array(spec["eigenvectors_desc_columns"])
        for col in (8, 9):
            comps = evecs[:, col]
            order = np.argsort(-np.abs(comps))
            top = ", ".join(f"{NAMES[k]}:{comps[k]:+.3f}" for k in order[:4])
            print(f"  evec {col+1} (ratio {spec['ratios_to_max'][col]:.3e}) {top}")
        print("  sloppy participation", dict(zip(NAMES, np.round(spec["sloppy_participation"], 3))))

    g_full = claim_log_ratio(th0, kind="full")
    b_early = burden_at(th0, T_CLAIM_EARLY)
    b_late = burden_at(th0, T_CLAIM_LATE)
    k_true = kappa_of(th0)

    # Good reduction: quasi-steady stress, same c. Evaluated on both maps.
    th_qss = qss_theta(th0)
    g_qss = claim_log_ratio(th_qss, kind="qss")
    chi_qss = {amap: cost_against_reference(th_qss, amap, mu[amap], kind="qss") for amap in ("B", "L")}
    print("claim g_full", g_full, "g_qss", g_qss, "chi_qss", chi_qss)
    print("gates", {k: (v["gated"], v["fraction_in_sloppy_subspace"]) for k, v in gates.items()})

    # Good reduction, second writing: also drop c when the map cannot see A.
    # Burden trajectory of the QSS model does not depend on c, so g is identical.
    g_drop_c = g_qss

    # Bad reduction at the generating parameter, and at a canyon point where c is ×20.
    th_bad = th0.copy()
    g_bad = claim_log_ratio(th_bad, kind="bad")
    chi_bad = {amap: cost_against_reference(th_bad, amap, mu[amap], kind="bad") for amap in ("B", "L")}
    th_c20 = th0.copy()
    th_c20[NAMES.index("c")] *= 20.0
    g_c20_full = claim_log_ratio(th_c20, kind="full")
    chi_c20 = {amap: cost_against_reference(th_c20, amap, mu[amap], kind="full") for amap in ("B", "L")}
    g_c20_bad = claim_log_ratio(th_c20, kind="bad")
    chi_c20_bad = {amap: cost_against_reference(th_c20, amap, mu[amap], kind="bad") for amap in ("B", "L")}

    factors = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    scans = {
        amap: scan_parameter(th0, amap, mu[amap], "c", factors) for amap in ("B", "L")
    }

    print("running eigenvector flow on map L (stress scale is the expected canyon)")
    flow_L = eigenvector_flow(th0, "L", mu["L"])
    print("  steps", len(flow_L) - 1, "final chi2", flow_L[-1]["chi2"])
    print("  final ratios", dict(zip(NAMES, np.round(flow_L[-1]["parameter_ratio_to_start"], 3))))

    print("running eigenvector flow on map B")
    flow_B = eigenvector_flow(th0, "B", mu["B"])
    print("  steps", len(flow_B) - 1, "final chi2", flow_B[-1]["chi2"])
    print("  final ratios", dict(zip(NAMES, np.round(flow_B[-1]["parameter_ratio_to_start"], 3))))

    # Hessian check. Noise-free data: 0.5 H should match Fisher.
    # Noisy data: one draw, same Jacobian point, residual nonzero.
    print("hessian checks")
    half_nf, evals_nf = observed_hessian(th0, "B", mu["B"])
    fim_B = fishers["B"]
    evals_f = np.sort(np.linalg.eigvalsh(fim_B))[::-1]
    rel = np.abs(evals_nf - evals_f) / np.maximum(np.abs(evals_f), 1e-8)
    noise = RNG.normal(0.0, SIGMA["B"], size=mu["B"].shape)
    # Four replicates: average of four noisy means is not required.
    # One series with variance sigma^2/n_rep is the sufficient statistic.
    # We put the noise on the mean at the replicate-averaged scale.
    mu_noisy = mu["B"] + noise / np.sqrt(N_REP["B"])
    half_ny, evals_ny = observed_hessian(th0, "B", mu_noisy)
    rel_ny = np.abs(evals_ny - evals_f) / np.maximum(np.abs(evals_f), 1e-8)

    # State samples for the table.
    dense_t = np.linspace(0.0, T_END, 401)
    Y = trajectory(th0, dense_t, kind="full")
    Yq = trajectory(th_qss, dense_t, kind="qss")
    Yb = trajectory(th_bad, dense_t, kind="bad")
    Yb20 = trajectory(th_c20, dense_t, kind="bad")

    def series_at(Ystate, times, has_S=True):
        # interpolate linearly on the dense grid
        rows = []
        for t in times:
            idx = int(np.argmin(np.abs(dense_t - t)))
            P, Q, A = Ystate[0, idx], Ystate[1, idx], Ystate[2, idx]
            row = {"t": float(t), "P": float(P), "Q": float(Q), "A": float(A), "B": float(P + Q)}
            if has_S and Ystate.shape[0] == 4:
                row["S"] = float(Ystate[3, idx])
            rows.append(row)
        return rows

    sample_times = [0.0, 11.0, 40.0]

    results = {
        "seed": SEED,
        "parameter_names": NAMES,
        "theta0": th0.tolist(),
        "y0": Y0.tolist(),
        "t_obs": T_OBS.tolist(),
        "sigma": SIGMA,
        "n_rep": N_REP,
        "cuts": {"numerical": NUM_CUT, "sloppy": SLOPPY_CUT},
        "kappa": k_true,
        "claim": {
            "definition": "g = log B(40) - log B(11), B = P + Q",
            "g_full": g_full,
            "B_11": b_early,
            "B_40": b_late,
            "ratio_B40_over_B11": b_late / b_early,
            "gates": gates,
        },
        "spectra": specs,
        "qss_reduction": {
            "description": "S = (gp/gs)*B and h*S = kappa*B with kappa = h*gp/gs. c retained in the three-state QSS model.",
            "theta_qss_names": ["r", "a", "aq", "K", "q", "f", "c", "kappa"],
            "theta_qss": th_qss.tolist(),
            "g_qss": g_qss,
            "abs_log_claim_error": abs(g_qss - g_full),
            "rel_claim_error": rel_err(g_qss, g_full),
            "chi2_vs_full_mean": chi_qss,
            "preserves_sign": bool(np.sign(g_qss) == np.sign(g_full)),
        },
        "bad_reduction": {
            "description": "dP/dt gains + rho*c*A with rho=0.35 fixed. c is then narrated as rebound biology.",
            "rho": BAD_RHO,
            "g_at_theta0": g_bad,
            "abs_log_claim_error_at_theta0": abs(g_bad - g_full),
            "chi2_at_theta0": chi_bad,
            "c_factor_20": {
                "g_full_model": g_c20_full,
                "chi2_full_model": chi_c20,
                "g_bad_model": g_c20_bad,
                "chi2_bad_model": chi_c20_bad,
                "abs_log_claim_error_bad": abs(g_c20_bad - g_full),
            },
        },
        "c_scans": scans,
        "flow_L": flow_L,
        "flow_B": flow_B,
        "hessian_map_B": {
            "noise_free_half_hessian_eigenvalues_desc": evals_nf.tolist(),
            "fisher_eigenvalues_desc": evals_f.tolist(),
            "max_rel_abs_diff_noise_free": float(np.max(rel)),
            "noisy_half_hessian_eigenvalues_desc": evals_ny.tolist(),
            "max_rel_abs_diff_noisy": float(np.max(rel_ny)),
            "noise_note": "One draw, seed 20260921, added to map B at scale sigma/sqrt(n_rep).",
        },
        "samples_full": series_at(Y, sample_times, True),
        "samples_qss": series_at(Yq, sample_times, False),
        "samples_bad": series_at(Yb, sample_times, True),
    }

    # Figures
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for amap, marker in (("B", "o"), ("L", "s")):
        ev = np.array(specs[amap]["eigenvalues_desc"])
        ax.semilogy(np.arange(1, 11), np.maximum(ev, 1e-16), marker=marker, label=f"map {amap}")
    ax.axhline(specs["B"]["lambda_max"] * SLOPPY_CUT, color="0.4", ls="--", lw=0.8, label="map B sloppy cut")
    ax.set_xlabel("Eigenvalue index (stiff to sloppy)")
    ax.set_ylabel("Fisher eigenvalue")
    ax.set_xticks(range(1, 11))
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fim_eigenspectra.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0), sharey=True)
    for ax, amap in zip(axes, ("B", "L")):
        part = np.array(specs[amap]["sloppy_participation"])
        ax.barh(NAMES[::-1], part[::-1], color="#3c3c3c")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Sloppy-subspace weight")
        ax.set_title(f"Map {amap}")
    fig.tight_layout()
    fig.savefig(FIG / "sloppy_participation.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for label, series, ls in (
        ("full", scans["B"], "-"),
        ("map L", scans["L"], "--"),
    ):
        fac = [row["factor"] for row in series]
        chi = [max(row["chi2"], 1e-18) for row in series]
        ax.loglog(fac, chi, ls=ls, marker="o", label=label if label != "full" else "map B")
    ax.set_xlabel("Factor applied to c, other rates fixed")
    ax.set_ylabel("Noise-free chi-square against the full mean")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "clearance_canyon.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=False)
    for ax, flow, title in (
        (axes[0], flow_B, "Map B flow"),
        (axes[1], flow_L, "Map L flow"),
    ):
        steps = [row["step"] for row in flow]
        for j, name in enumerate(NAMES):
            series = [row["parameter_ratio_to_start"][j] for row in flow]
            ax.plot(steps, series, label=name, lw=1.2)
        ax.set_yscale("log")
        ax.set_ylabel("Parameter / start")
        ax.set_title(title)
        ax.axhline(1.0, color="0.5", lw=0.6)
    axes[1].set_xlabel("Flow step")
    axes[0].legend(ncol=5, frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "mbam_flow.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    Bfull = Y[0] + Y[1]
    Bq = Yq[0] + Yq[1]
    Bbad = Yb[0] + Yb[1]
    Bbad20 = Yb20[0] + Yb20[1]
    ax.plot(dense_t, Bfull, color="black", label="full model")
    ax.plot(dense_t, Bq, color="#1f4e79", ls="--", label="QSS stress reduction")
    ax.plot(dense_t, Bbad, color="#8c2f39", ls=":", label="bad rebound at generating c")
    ax.plot(dense_t, Bbad20, color="#8c2f39", ls="-.", label="bad rebound at 20× c")
    ax.axvline(T_CLAIM_EARLY, color="0.6", lw=0.6)
    ax.axvline(T_CLAIM_LATE, color="0.6", lw=0.6)
    ax.set_xlabel("Model time")
    ax.set_ylabel("Viable burden B = P + Q")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "burden_reductions.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    labels = ["full", "QSS\nreduction", "bad model\nat c", "bad model\nat 20× c"]
    vals = [g_full, g_qss, g_bad, g_c20_bad]
    ax.bar(labels, vals, color=["#222222", "#1f4e79", "#8c2f39", "#c47b84"])
    ax.axhline(0.0, color="0.4", lw=0.6)
    ax.set_ylabel("g = log B(40) - log B(11)")
    fig.tight_layout()
    fig.savefig(FIG / "claim_values.png", dpi=160)
    plt.close(fig)

    def flow_endpoint(flow):
        row = flow[-1]
        th = th0 * np.array(row["parameter_ratio_to_start"])
        g_end = claim_log_ratio(th, kind="full")
        return {
            "steps": int(row["step"]),
            "chi2": row["chi2"],
            "g": g_end,
            "abs_g_change": abs(g_end - g_full),
            "kappa": row["kappa"],
            "parameter_ratio_to_start": row["parameter_ratio_to_start"],
        }

    def rms_from_chi2(chi, amap, n_obs):
        weight = N_REP[amap] / (SIGMA[amap] ** 2)
        return float(np.sqrt(chi / weight / n_obs))

    sloppy_L = np.array(specs["L"]["eigenvectors_desc_columns"])
    # columns are stiff to sloppy; last two are the stress plane
    idx = {name: i for i, name in enumerate(NAMES)}
    dlog_kappa = []
    for col in (8, 9):
        v = sloppy_L[:, col]
        dlog = float(v[idx["h"]] + v[idx["gp"]] - v[idx["gs"]])
        dlog_kappa.append({"column_stiff_to_sloppy": col + 1, "d_log_kappa": dlog})

    bad_sweep = []
    for fac in (0.05, 0.2, 1.0, 5.0, 20.0):
        th = th0.copy()
        th[idx["c"]] *= fac
        bad_sweep.append(
            {
                "factor": fac,
                "g_full": claim_log_ratio(th, kind="full"),
                "chi2_full_B": cost_against_reference(th, "B", mu["B"], kind="full"),
                "g_bad": claim_log_ratio(th, kind="bad"),
                "chi2_bad_B": cost_against_reference(th, "B", mu["B"], kind="bad"),
            }
        )

    fish5 = np.array(results["hessian_map_B"]["fisher_eigenvalues_desc"][:5])
    hnf5 = np.array(results["hessian_map_B"]["noise_free_half_hessian_eigenvalues_desc"][:5])
    hny5 = np.array(results["hessian_map_B"]["noisy_half_hessian_eigenvalues_desc"][:5])
    results["hessian_map_B"]["rel_abs_diff_leading5_noise_free"] = (
        np.abs(hnf5 - fish5) / np.abs(fish5)
    ).tolist()
    results["hessian_map_B"]["rel_abs_diff_leading5_noisy"] = (
        np.abs(hny5 - fish5) / np.abs(fish5)
    ).tolist()

    def _only(theta, row, name):
        out_th = np.array(theta, dtype=float)
        j = NAMES.index(name)
        out_th[j] = theta[j] * row["parameter_ratio_to_start"][j]
        return out_th

    results["summary"] = {
        "flow_B_endpoint": flow_endpoint(flow_B),
        "flow_L_endpoint": flow_endpoint(flow_L),
        "qss_rms_log": {"B": rms_from_chi2(chi_qss["B"], "B", T_OBS.size), "L": rms_from_chi2(chi_qss["L"], "L", 3 * T_OBS.size)},
        "map_L_sloppy_d_log_kappa": dlog_kappa,
        "bad_c_sweep": bad_sweep,
        "only_a_at_flow_B_endpoint_chi2": cost_against_reference(
            _only(th0, flow_B[-1], "a"), "B", mu["B"]
        ),
    }

    out = ROOT / "results.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("wrote", out)
    print("summary", json.dumps(results["summary"], indent=2)[:2500])


if __name__ == "__main__":
    main()
