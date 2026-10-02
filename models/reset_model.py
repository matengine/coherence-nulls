"""
Phenomenological two-region model in which the fast rhythms never interact.

Each region has its own gamma oscillator with its own frequency jitter. The two share a theta rhythm that
partially resets each oscillator at every trough, and nothing else; `coupling_c` is an optional positive
control that transmits region A's gamma into region B's field. Fields are continuous, not trials, so they
go straight into coherence.spectra.

Three parameters carry the argument and are all measurable in data:

  reset_r               how completely a trough resets the local phase; weak resets give no coherence
  reset_time_jitter_ms  independent jitter of the two regions' reset times; 2 ms costs about half the
                        coherence, 4 ms all of it
  theta_local_amp       an independent theta-band component added to each field only, which lowers the
                        measured theta coherence without touching the resets. Lowering theta coherence the
                        other way, by letting the two regions' theta phases wander apart
                        (theta_local_sd_deg), also destroys the gamma coherence, because the resets stop
                        coinciding. The two routes make different predictions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from models.noise import pink_noise, ou_process
from coherence.circular_statistics import wrap, troughs_from_phase


@dataclass
class ModelParams:
    fs: float = 1250.0
    duration_s: float = 500.0

    # shared slow rhythm
    theta_mean_hz: float = 8.8
    theta_drift_sd_hz: float = 0.4          # 0.1 Hz with tau 5 s is the "stable theta" regime
    theta_drift_tau_s: float = 2.0
    theta_amp: float = 1.6
    theta_offset_deg: float = 60.0          # fixed anatomical phase offset of region B

    # route (i): each region's slow phase wanders, which moves its reset times
    theta_local_sd_deg: float = 0.0
    theta_local_tau_s: float = 1.0

    # route (ii): an independent slow-band component in each FIELD, which does not move the resets
    theta_local_amp: float = 0.0
    theta_local_drift_sd_hz: float = 0.1    # -1 means the same drift as the shared rhythm
    theta_local_drift_tau_s: float = 5.0
    theta_local_f_offset_hz: float = 0.0

    # non-sinusoidal slow waveform: [(harmonic_number, amplitude, phase_rad), ...]
    theta_harmonics: Sequence[tuple] = field(default_factory=tuple)

    # local fast oscillators
    gamma_f0: Sequence[float] = (72.0, 78.0)
    gamma_amp: float = 0.70
    gamma_jitter_hz: float = 5.0
    gamma_jitter_tau_s: float = 0.03
    gamma_env_depth: float = 0.8            # slow-rhythm amplitude modulation of the fast rhythm

    # the reset
    reset_r: float = 0.32                   # fraction of the phase error removed at each trough
    reset_time_jitter_ms: float = 0.0       # independent jitter of each region's reset times

    noise_amp: float = 0.6

    # positive control: direct transmission of A's fast signal into B's field
    coupling_c: float = 0.0
    coupling_delay_ms: float = 5.0

    seed: int = 0


def _reset_oscillator(n, fs, f0, jitter_hz, jitter_tau, reset_idx, r, rng):
    """A free-running oscillator whose phase is pulled a fraction r toward zero at each reset index."""
    finst = f0 + ou_process(n, fs, jitter_hz, jitter_tau, rng)
    phi_free = np.cumsum(2 * np.pi * finst / fs) + rng.uniform(0, 2 * np.pi)
    J = 0.0
    jumps = np.zeros(len(reset_idx))
    for j, k in enumerate(reset_idx):
        cur = phi_free[k] + J
        J += r * wrap(0.0 - cur)
        jumps[j] = J
    step = np.zeros(n)
    if len(reset_idx):
        step[reset_idx] = np.diff(np.concatenate([[0.0], jumps]))
    return phi_free + np.cumsum(step)


def simulate(p: ModelParams):
    """Run the model. Returns a dict with 'lfp_A', 'lfp_B' and the latent variables behind them.

    Keys: t, fs, f_theta, theta_A, theta_B, troughs_A/B, gamma_phase_A/B, gamma_A/B, lfp_A, lfp_B, and
    theta_local_A/B when an independent slow-band component is present. The troughs are the true reset
    times, which is what a conditional null needs and what a recording does not supply.
    """
    rng = np.random.default_rng(p.seed)
    fs = p.fs
    n = int(round(p.duration_s * fs))
    t = np.arange(n) / fs

    f_theta = np.clip(p.theta_mean_hz + ou_process(n, fs, p.theta_drift_sd_hz, p.theta_drift_tau_s, rng), 5.0, 11.0)
    theta_shared = np.cumsum(2 * np.pi * f_theta / fs) + rng.uniform(0, 2 * np.pi)

    sd_loc = np.radians(p.theta_local_sd_deg)
    psi_A = ou_process(n, fs, sd_loc, p.theta_local_tau_s, rng) if sd_loc > 0 else np.zeros(n)
    psi_B = ou_process(n, fs, sd_loc, p.theta_local_tau_s, rng) if sd_loc > 0 else np.zeros(n)
    theta_A = theta_shared + psi_A
    theta_B = theta_shared + psi_B + np.radians(p.theta_offset_deg)

    def theta_wave(th):
        w = np.cos(th)
        for (k, amp, ph) in p.theta_harmonics:
            w = w + amp * np.cos(k * th + ph)
        return w

    out = {"t": t, "fs": fs, "f_theta": f_theta, "theta_A": theta_A, "theta_B": theta_B}
    gam = {}
    for name, th, f0 in (("A", theta_A, p.gamma_f0[0]), ("B", theta_B, p.gamma_f0[1])):
        troughs = troughs_from_phase(th)
        if p.reset_time_jitter_ms > 0:
            jit = np.round(rng.normal(0, p.reset_time_jitter_ms * fs / 1000.0, len(troughs))).astype(int)
            troughs = np.unique(np.clip(troughs + jit, 0, n - 1))
        phi = _reset_oscillator(n, fs, f0, p.gamma_jitter_hz, p.gamma_jitter_tau_s, troughs, p.reset_r, rng)
        env = 1.0 + p.gamma_env_depth * np.cos(th - np.pi)
        g = env * np.cos(phi)
        gam[name] = g
        out["troughs_" + name] = troughs
        out["gamma_phase_" + name] = phi
        out["gamma_" + name] = g

    lfp_A = p.theta_amp * theta_wave(theta_A) + p.gamma_amp * gam["A"] + p.noise_amp * pink_noise(n, rng)
    lfp_B = p.theta_amp * theta_wave(theta_B) + p.gamma_amp * gam["B"] + p.noise_amp * pink_noise(n, rng)

    if p.theta_local_amp > 0:
        sd_l = p.theta_local_drift_sd_hz if p.theta_local_drift_sd_hz >= 0 else p.theta_drift_sd_hz
        tau_l = p.theta_local_drift_tau_s if p.theta_local_drift_tau_s > 0 else p.theta_drift_tau_s
        for name in ("A", "B"):
            off = -p.theta_local_f_offset_hz if name == "A" else p.theta_local_f_offset_hz
            f_loc = np.clip(p.theta_mean_hz + off + ou_process(n, fs, sd_l, tau_l, rng), 5.0, 11.0)
            th_loc = np.cumsum(2 * np.pi * f_loc / fs) + rng.uniform(0, 2 * np.pi)
            out["theta_local_" + name] = th_loc
            if name == "A":
                lfp_A = lfp_A + p.theta_local_amp * p.theta_amp * theta_wave(th_loc)
            else:
                lfp_B = lfp_B + p.theta_local_amp * p.theta_amp * theta_wave(th_loc)

    if p.coupling_c > 0:
        d = int(round(p.coupling_delay_ms * fs / 1000.0))
        transmitted = np.zeros(n)
        transmitted[d:] = gam["A"][: n - d]
        lfp_B = lfp_B + p.coupling_c * p.gamma_amp * transmitted

    out["lfp_A"] = lfp_A
    out["lfp_B"] = lfp_B
    return out


# Speed-linked operating points. Only the reset reliability and the slow-rhythm frequency are tied to
# running speed, both of which rise with speed in the hippocampus; every other parameter is held fixed, so
# any speed dependence of the resulting coherence is attributable to those two alone.
SPEEDS_CM_S = (2.6, 8.7, 24.0, 46.0)
SPEED_RESET = (0.10, 0.20, 0.32, 0.45)
SPEED_THETA_HZ = (7.4, 7.8, 8.3, 8.8)


def speed_params(i, **kw):
    """ModelParams for speed bin i (0 slowest, 3 fastest), with any field overridden by keyword."""
    p = ModelParams(reset_r=SPEED_RESET[i], theta_mean_hz=SPEED_THETA_HZ[i], seed=100 + i)
    for k, v in kw.items():
        setattr(p, k, v)
    return p
