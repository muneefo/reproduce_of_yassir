"""Section IV parameters; every physical quantity is in SI units."""
from dataclasses import asdict, dataclass, replace
import numpy as np


@dataclass(frozen=True)
class Config:
    density_m2: float = 70 / 1e6
    carrier_hz: float = 28e9
    bandwidth_hz: float = 200e6
    beams: int = 8
    beamwidth_rad: float = np.deg2rad(18)
    gain: float = 10 ** (20 / 10)
    alpha_l: float = 2
    alpha_n: float = 4
    intercept_l: float = 10 ** (-61.4 / 10)
    intercept_n: float = 10 ** (-72 / 10)
    m_l: float = 3
    m_n: float = 2
    rcs_mono_m2: float = 1
    rcs_bi_m2: float = .7
    blockage_m_inv: float = .0149
    average_power_w: float = 1
    gamma: float = .7
    temperature_k: float = 300
    residual_si: float = 1e-12
    speed_m_s: float = 299792458
    boltzmann_j_k: float = 1.380649e-23
    area_m2: float = 25e6
    seed: int = 20261010
    # Choice, not specified in the paper: budget interpreted per beam.
    power_scope: str = "per_beam"

    def __post_init__(self):
        if not 0 <= self.gamma <= 1 or self.density_m2 <= 0:
            raise ValueError("Invalid energy split or density")
        if self.power_scope not in ("per_beam", "per_bs"):
            raise ValueError("Unknown power scope")
        if self.beams * 2 * self.beamwidth_rad > 2 * np.pi:
            raise ValueError("Chosen equally spaced cosine supports overlap")
        if self.slot_s <= self.pulse_s:
            raise ValueError("Echo-reception interval must be positive")

    @property
    def wavelength_m(self):
        return self.speed_m_s / self.carrier_hz

    @property
    def range_m(self):
        return 3 / np.sqrt(np.pi * self.density_m2)

    @property
    def pulse_s(self):
        return 1 / self.bandwidth_hz

    @property
    def slot_s(self):
        return 2 * self.range_m / self.speed_m_s

    @property
    def energy_j(self):
        divisor = self.beams if self.power_scope == "per_bs" else 1
        return self.average_power_w * self.slot_s / divisor

    @property
    def sensing_w(self):
        return self.gamma * self.energy_j / self.pulse_s

    @property
    def communication_w(self):
        return (1 - self.gamma) * self.energy_j / (self.slot_s - self.pulse_s)

    @property
    def noise_w(self):
        return self.boltzmann_j_k * self.temperature_k * self.bandwidth_hz

    @property
    def radar_constant(self):
        return self.sensing_w * self.gain**2 * self.wavelength_m**2 / (4*np.pi)**3

    def changed(self, **kwargs):
        return replace(self, **kwargs)

    def to_dict(self):
        d = asdict(self)
        d.update({k: getattr(self, k) for k in
                  ("wavelength_m", "range_m", "pulse_s", "slot_s", "energy_j",
                   "sensing_w", "communication_w", "noise_w")})
        return d
