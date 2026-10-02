"""
Stochastic building blocks shared by the models.
"""
import numpy as np
from scipy import signal


def pink_noise(n, rng):
    """Unit-variance 1/f noise of length n, built by randomizing the phase of a 1/sqrt(f) amplitude spectrum."""
    f = np.fft.rfftfreq(n)
    amp = np.zeros_like(f)
    amp[1:] = 1.0 / np.sqrt(f[1:])
    spec = amp * np.exp(1j * rng.uniform(0, 2 * np.pi, len(f)))
    x = np.fft.irfft(spec, n)
    return (x - x.mean()) / x.std()


def ou_process(n, fs, sd, tau, rng):
    """Ornstein-Uhlenbeck process with stationary SD `sd` and time constant `tau` in seconds.

    Used for slow drift of a rhythm's frequency. A perfectly periodic rhythm is a trap: it locks two
    independent oscillators into exact harmonic relationships and produces spurious near-unity coherence at
    the harmonics, so the frequency must be allowed to vary.
    """
    a = np.exp(-1.0 / (fs * tau))
    e = rng.normal(0.0, sd * np.sqrt(1.0 - a ** 2), n)
    return signal.lfilter([1.0], [1.0, -a], e, zi=[rng.normal(0.0, sd)])[0]
