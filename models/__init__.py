"""
The two generative models the analyses are validated against.

  reset_model    two gamma oscillators that never interact, sharing only a theta rhythm that partially
                 resets each. Reproduces the reported cross-regional gamma coherence with no gamma
                 interaction, and shows what that costs in required timing precision.

  ping_network   two pyramidal-interneuron gamma networks of leaky integrate-and-fire neurons, with a
                 field proxy from the recurrent synaptic currents onto the pyramidal cells. A biophysical
                 generator does not behave like a free-running oscillator: smooth shared modulation that
                 leaves it running gives no gamma coherence, and a shared reset aligns the two rhythms
                 only for their own phase memory, one or two cycles.

  noise          pink noise and the Ornstein-Uhlenbeck process used for slow drift.

Use these as the harness for any null meant for real data: run it where the answer is known first.
"""
from .noise import pink_noise, ou_process
from .reset_model import (ModelParams, simulate, speed_params, SPEEDS_CM_S, SPEED_RESET, SPEED_THETA_HZ)
from .ping_network import PingParams, simulate_ping

__all__ = ["pink_noise", "ou_process", "ModelParams", "simulate", "speed_params",
           "SPEEDS_CM_S", "SPEED_RESET", "SPEED_THETA_HZ", "PingParams", "simulate_ping"]
