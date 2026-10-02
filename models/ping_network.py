"""
Two independent PING (pyramidal-interneuron gamma) networks of leaky integrate-and-fire neurons that
share only an external theta-paced drive, with a synaptic-current field proxy.

Design choices (with sources):
  * LIF membranes and fast synapses in the Brunel & Wang (2003) regime: tau_m = 20 ms (E) / 10 ms (I),
    AMPA tau = 2 ms, GABA_A tau = 4 ms, 0.75 ms synaptic delay, refractory 2 ms (E) / 1 ms (I).
  * Current-based exponential synapses; the network sits in the fluctuation-driven regime with
    independent external Poisson input to every neuron, so that each region's gamma phase is set
    by its own noise unless the shared drive resets it.
  * Field proxy per region: the sum over pyramidal cells of the absolute AMPA and GABA synaptic
    currents (Mazzoni et al., 2015), sampled at 1250 Hz. Two versions are returned: 'recurrent'
    (local E->E AMPA and I->E GABA only; the default, so that the field reflects the network's own
    activity) and 'full' (recurrent plus the AMPA current of the external drive itself, which carries
    the shared drive's waveform into the field with no network involvement).
  * Shared drive modes: 'tonic' (constant external rate), 'smooth' (rate multiplied by 1 + m cos theta),
    'pulsed' (same mean rate, with the modulated part delivered as a brief Gaussian barrage once per
    theta cycle at the theta trough). Region B's barrage times can be jittered relative to region A's
    to test the timing requirement.
  * No synapses of any kind between the two regions; an optional A->B transmission term for validation.

The simulation loop is plain numpy with dt = 0.1 ms.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy import signal


@dataclass
class PingParams:
    # populations
    NE: int = 400
    NI: int = 100
    # membrane (mV, ms)
    tau_mE: float = 20.0
    tau_mI: float = 10.0
    EL: float = -70.0
    Vth: float = -52.0
    Vreset: float = -60.0
    tref_E: float = 2.0
    tref_I: float = 1.0
    # Synapses (current-based exponential; w = total depolarizing charge per spike in mV).
    # With these defaults and tonic drive the network settles into the
    # sparsely synchronized regime, a 68 Hz population rhythm with a quality factor near 2.6 and a phase
    # memory of one to two cycles, while each pyramidal cell fires on about one gamma cycle in seven.
    # Two nearby settings leave that regime and are worth knowing about: lowering the input noise, or
    # lengthening the synaptic delay past a millisecond, pushes the network into full synchrony, a 30-50 Hz
    # rhythm on which a third of the pyramidal cells fire every cycle and which carries a harmonic comb of
    # its own; and the sharpness of a sparsely synchronized rhythm is limited by finite-size fluctuations,
    # so larger networks with the same mean coupling give sharper rhythms and longer alignment after a reset.
    tau_A: float = 2.0
    tau_G: float = 4.0
    delay_ms: float = 0.75
    p_EE: float = 0.10
    p_EI: float = 0.20
    p_IE: float = 0.20
    p_II: float = 0.20
    w_EE: float = 0.20
    w_EI: float = 0.60
    w_IE: float = 2.00
    w_II: float = 0.80
    # external drive: independent Poisson input to every neuron (rate in kHz, weight in mV)
    ext_rate_E: float = 3.0
    ext_rate_I: float = 3.0
    w_ext: float = 0.5
    # shared theta drive
    mode: str = "pulsed"        # 'tonic' | 'smooth' | 'pulsed'
    theta_mean_hz: float = 8.8            # as the phenomenological model's 46 cm/s operating point
    theta_drift_sd_hz: float = 0.4        # same drifting-theta regime as the phenomenological model
    theta_drift_tau_s: float = 2.0
    theta_offset_deg: float = 60.0
    mod_depth: float = 0.6      # fraction of the mean external rate that is theta-modulated
    drive_target: str = "all"   # which cells receive the theta modulation: 'all' or 'E' (interneurons then tonic)
    # 'mixed' mode: shared smooth modulation (mod_depth) + an independent local theta modulation per network
    # (local_depth; own phase, drift local_drift_sd_hz / local_drift_tau_s) + an optional volley of per-cell size
    # volley_frac (fraction of the mean rate, Gaussian SD volley_sd_ms) at the shared theta trough, delivered to a
    # random fraction volley_subset of the excitatory cells only
    local_depth: float = 0.0        # private theta modulating the network's DRIVE (depth of the rate modulation)
    local_field_amp: float = 0.0    # private theta added to the FIELD only (RMS relative to the field's own 5-11 Hz RMS),
                                    # the counterpart of the phenomenological model's independent theta-band component
    local_drift_sd_hz: float = 0.1
    local_drift_tau_s: float = 5.0
    volley_frac: float = 0.0
    volley_sd_ms: float = 3.0
    volley_subset: float = 1.0
    volley_target: str = "E"    # 'E': extra excitation to a subset of E cells; 'I': extra excitation to a subset of the
                                # interneurons (a synchronous IPSP on the pyramidal cells); 'Epause': a transient
                                # reduction of the E cells' drive (rate clipped at zero)
    pulse_sd_ms: float = 3.0    # width of the barrage (Gaussian SD) in the pulsed mode
    barrage_jitter_ms: float = 0.0   # jitter of region B's barrage times relative to region A's
    # optional direct transmission A -> B (fraction of region A's E-population spikes delivered to B's E and I cells)
    coupling_c: float = 0.0
    coupling_delay_ms: float = 5.0
    # field proxy returned as lfp_A / lfp_B: 'recurrent' (default) or 'full' (see module docstring)
    field: str = "recurrent"
    # simulation
    dt: float = 0.1             # ms
    duration_s: float = 20.0
    fs_out: float = 1250.0
    seed: int = 0
    record_raster_s: float = 1.0


def _ou(n, fs, sd, tau, rng):
    a = np.exp(-1.0 / (fs * tau))
    e = rng.normal(0.0, sd * np.sqrt(1.0 - a ** 2), n)
    return signal.lfilter([1.0], [1.0, -a], e, zi=[rng.normal(0.0, sd)])[0]


def _troughs(theta):
    w = np.angle(np.exp(1j * (theta - np.pi)))
    return np.where((w[:-1] < 0) & (w[1:] >= 0))[0] + 1


def simulate_ping(p: PingParams):
    rng = np.random.default_rng(p.seed)
    dt = p.dt
    nsteps = int(round(p.duration_s * 1000.0 / dt))
    N = p.NE + p.NI                     # per region
    NT = 2 * N                          # both regions in one state vector
    isE = np.zeros(NT, bool); isE[:p.NE] = True; isE[N:N + p.NE] = True
    region = np.zeros(NT, int); region[N:] = 1

    # ---- connectivity (block diagonal; weights already divided by tau so that s integrates to w)
    def block(nrows, ncols, prob, w):
        M = (rng.random((nrows, ncols)) < prob).astype(np.float32) * np.float32(w)
        return M
    W = np.zeros((NT, NT), np.float32)   # W[post, pre]
    for r in range(2):
        o = r * N
        E = slice(o, o + p.NE); I = slice(o + p.NE, o + N)
        W[E, E] = block(p.NE, p.NE, p.p_EE, p.w_EE); np.fill_diagonal(W[E, E], 0)
        W[I, E] = block(p.NI, p.NE, p.p_EI, p.w_EI)
        W[E, I] = block(p.NE, p.NI, p.p_IE, -p.w_IE)
        W[I, I] = block(p.NI, p.NI, p.p_II, -p.w_II); np.fill_diagonal(W[I, I], 0)
    if p.coupling_c > 0:   # A's E cells -> B's E and I cells, same statistics as local E projections, scaled by c
        EA = slice(0, p.NE); EB = slice(N, N + p.NE); IB = slice(N + p.NE, 2 * N)
        W[EB, EA] = block(p.NE, p.NE, p.p_EE, p.w_EE * p.coupling_c)
        W[IB, EA] = block(p.NI, p.NE, p.p_EI, p.w_EI * p.coupling_c)
    Wpos = np.where(W > 0, W, 0).astype(np.float32)
    Wneg = np.where(W < 0, -W, 0).astype(np.float32)

    # ---- shared theta drive
    fs_sim = 1000.0 / dt
    f_theta = np.clip(p.theta_mean_hz + _ou(nsteps, fs_sim, p.theta_drift_sd_hz, p.theta_drift_tau_s, rng), 5.0, 11.0)
    theta = np.cumsum(2 * np.pi * f_theta / fs_sim) + rng.uniform(0, 2 * np.pi)
    theta_B = theta + np.radians(p.theta_offset_deg)
    rate_mult = np.ones((2, nsteps), np.float32)
    barrage_A = _troughs(theta); barrage_B = _troughs(theta_B)
    if p.barrage_jitter_ms > 0:
        barrage_B = np.unique(np.clip(barrage_B + np.round(rng.normal(0, p.barrage_jitter_ms / dt, len(barrage_B))).astype(int), 0, nsteps - 1))
    if p.mode == "smooth":
        rate_mult[0] = 1.0 + p.mod_depth * np.cos(theta - np.pi)
        rate_mult[1] = 1.0 + p.mod_depth * np.cos(theta_B - np.pi)
    elif p.mode == "pulsed":
        sd_steps = p.pulse_sd_ms / dt
        kern_half = int(4 * sd_steps)
        kk = np.arange(-kern_half, kern_half + 1)
        kern = np.exp(-0.5 * (kk / sd_steps) ** 2); kern /= kern.sum()
        for r, bt in enumerate((barrage_A, barrage_B)):
            imp = np.zeros(nsteps, np.float32); imp[bt] = 1.0
            pulses = np.convolve(imp, kern, mode="same")
            mean_cycle_steps = nsteps / max(len(bt), 1)
            rate_mult[r] = (1.0 - p.mod_depth) + p.mod_depth * pulses * mean_cycle_steps
    volley_mult = np.zeros((2, nsteps), np.float32)     # extra multiplier for the volley-receiving E cells ('mixed' mode)
    vol_mask = np.zeros(NT, np.float32)
    th_locs = []                                        # private theta phase per region (drive and/or field component)
    for r in range(2):
        f_loc = np.clip(p.theta_mean_hz + _ou(nsteps, fs_sim, p.local_drift_sd_hz, p.local_drift_tau_s, rng), 5.0, 11.0)
        th_locs.append(np.cumsum(2 * np.pi * f_loc / fs_sim) + rng.uniform(0, 2 * np.pi))
    if p.mode == "mixed":
        for r, (th, bt) in enumerate(((theta, barrage_A), (theta_B, barrage_B))):
            th_loc = th_locs[r]
            rate_mult[r] = 1.0 + p.mod_depth * np.cos(th - np.pi) + p.local_depth * np.cos(th_loc - np.pi)
            if p.volley_frac > 0:
                sd_steps = p.volley_sd_ms / dt; kk = np.arange(-int(4 * sd_steps), int(4 * sd_steps) + 1)
                kern = np.exp(-0.5 * (kk / sd_steps) ** 2); kern /= kern.sum()
                imp = np.zeros(nsteps, np.float32); imp[bt] = 1.0
                volley_mult[r] = p.volley_frac * np.convolve(imp, kern, mode="same") * (nsteps / max(len(bt), 1))
        for r in range(2):
            o = r * N
            if p.volley_target == "I":
                idx = p.NE + rng.choice(p.NI, int(round(p.volley_subset * p.NI)), replace=False)
            else:
                idx = rng.choice(p.NE, int(round(p.volley_subset * p.NE)), replace=False)
            vol_mask[o + idx] = -1.0 if p.volley_target == "Epause" else 1.0
    # transmission delay ring buffer
    D = max(int(round(p.delay_ms / dt)), 1)
    Dc = max(int(round(p.coupling_delay_ms / dt)), 1)

    # ---- state
    V = p.EL + rng.uniform(0, 10, NT)
    sA = np.zeros(NT, np.float32); sG = np.zeros(NT, np.float32)   # synaptic currents (mV/ms)
    sA_ext = np.zeros(NT, np.float32)                               # the external-drive part of sA (for the proxy)
    ref = np.zeros(NT, np.float32)
    tau_m = np.where(isE, p.tau_mE, p.tau_mI).astype(np.float32)
    tref = np.where(isE, p.tref_E, p.tref_I).astype(np.float32)
    ext_rate = np.where(isE, p.ext_rate_E, p.ext_rate_I).astype(np.float32) * dt  # expected external spikes per step
    mod_mask = (isE if p.drive_target == "E" else np.ones(NT, bool)).astype(np.float32)   # cells whose rate is modulated
    decA = np.float32(np.exp(-dt / p.tau_A)); decG = np.float32(np.exp(-dt / p.tau_G))
    spike_buf = np.zeros((D, NT), bool)
    coup_buf = np.zeros((Dc, NT), bool) if p.coupling_c > 0 else None

    # ---- output
    dec = int(round(fs_sim / p.fs_out))
    nout = nsteps // dec
    lfp = np.zeros((2, 2, nout), np.float32)   # [kind (0 recurrent, 1 full), region]
    rates = np.zeros((2, 2))   # [region, E/I] mean rate
    n_rast = int(p.record_raster_s * 1000 / dt)
    rast_t, rast_i = [], []
    wext_over_tauA = np.float32(p.w_ext / p.tau_A)
    inv_tauA = np.float32(1.0 / p.tau_A); inv_tauG = np.float32(1.0 / p.tau_G)
    acc_lfp = np.zeros((2, 2), np.float64); acc_n = 0; out_idx = 0
    spk_count = np.zeros(NT)

    for t in range(nsteps):
        # delayed recurrent spikes arriving now
        arriving = spike_buf[t % D]
        if arriving.any():
            idx = np.flatnonzero(arriving)
            sA += Wpos[:, idx].sum(axis=1) * inv_tauA
            sG += Wneg[:, idx].sum(axis=1) * inv_tauG
        if coup_buf is not None:
            arr_c = coup_buf[t % Dc]
            if arr_c.any():
                idx = np.flatnonzero(arr_c)
                sA += Wpos[:, idx].sum(axis=1) * inv_tauA
        # external Poisson input (theta-modulated per region)
        lam = np.maximum(ext_rate * (1.0 + (rate_mult[region, t] - 1.0) * mod_mask + volley_mult[region, t] * vol_mask), 0.0)
        next_ = rng.poisson(lam).astype(np.float32)
        sA += next_ * wext_over_tauA
        sA_ext += next_ * wext_over_tauA
        # membrane update (current-based, Euler)
        dV = (p.EL - V) / tau_m + (sA - sG)
        V = np.where(ref > 0, V, V + dt * dV)
        ref = np.maximum(ref - dt, 0)
        spk = V >= p.Vth
        if spk.any():
            V[spk] = p.Vreset; ref[spk] = tref[spk]
            spk_count += spk
            if t < n_rast:
                ii = np.flatnonzero(spk); rast_t.append(np.full(len(ii), t * dt)); rast_i.append(ii)
        # store spikes for delayed delivery (local synapses)
        local_spk = spk.copy()
        if coup_buf is not None:
            coup_buf[(t + Dc) % Dc] = False
            cs = np.zeros(NT, bool); cs[:p.NE] = spk[:p.NE]     # only A's E cells cross over
            coup_buf[(t + Dc - 1) % Dc] = cs
        spike_buf[(t + D - 1) % D] = local_spk
        # synaptic decay
        sA *= decA; sG *= decG; sA_ext *= decA
        # field proxy: sum over E cells of |I_AMPA| + |I_GABA|; recurrent-only (row 0) and with the external current (row 1)
        gA = float(sG[:p.NE].sum()); gB = float(sG[N:N + p.NE].sum())
        aA = float(sA[:p.NE].sum()); aB = float(sA[N:N + p.NE].sum())
        eA = float(sA_ext[:p.NE].sum()); eB = float(sA_ext[N:N + p.NE].sum())
        acc_lfp[0, 0] += aA - eA + gA; acc_lfp[0, 1] += aB - eB + gB
        acc_lfp[1, 0] += aA + gA; acc_lfp[1, 1] += aB + gB
        acc_n += 1
        if acc_n == dec:
            if out_idx < nout:
                lfp[:, :, out_idx] = acc_lfp / dec
            out_idx += 1; acc_lfp[:] = 0; acc_n = 0

    if p.local_field_amp > 0:   # private theta-band component in the field only: RMS = local_field_amp x the field's 5-11 Hz RMS
        from scipy import signal as _sg
        bb, ab = _sg.butter(2, (5.0, 11.0), btype="bandpass", fs=p.fs_out)
        for r in range(2):
            loc = np.cos(th_locs[r][::dec][:nout] - np.pi).astype(np.float32)
            for kind in range(2):
                x = lfp[kind, r]
                rms_th = float(np.std(_sg.filtfilt(bb, ab, x - x.mean())))
                lfp[kind, r] = x + np.float32(p.local_field_amp * rms_th * np.sqrt(2.0)) * loc

    T = p.duration_s
    rates[0, 0] = spk_count[:p.NE].mean() / T; rates[0, 1] = spk_count[p.NE:N].mean() / T
    rates[1, 0] = spk_count[N:N + p.NE].mean() / T; rates[1, 1] = spk_count[N + p.NE:].mean() / T
    k = 1 if p.field == "full" else 0
    out = dict(lfp_A=lfp[k, 0].astype(float), lfp_B=lfp[k, 1].astype(float), fs=p.fs_out, rates=rates,
               lfp_A_rec=lfp[0, 0].astype(float), lfp_B_rec=lfp[0, 1].astype(float),
               lfp_A_full=lfp[1, 0].astype(float), lfp_B_full=lfp[1, 1].astype(float),
               theta=theta[::dec][:nout], theta_B=theta_B[::dec][:nout],
               barrage_A=barrage_A // dec, barrage_B=barrage_B // dec,
               raster_t=np.concatenate(rast_t) if rast_t else np.array([]),
               raster_i=np.concatenate(rast_i) if rast_i else np.array([]))
    return out
