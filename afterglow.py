from astropy import constants as const
import numpy as np

m_e = const.m_e.to("g").value  # type:ignore
m_p = const.m_p.to("g").value  # type:ignore
pi = np.pi
c = const.c.to("cm/s").value  # type:ignore
sigma_T = const.sigma_T.cgs.value  # type: ignore
q_e = const.e.esu.value  # type: ignore


class Afterglow:
    def __init__(self, E0, n, g0, eps_e, eps_b, p, d_l):
        self.E0 = E0
        self.n = n
        self.g0 = g0
        self.eps_e = eps_e
        self.eps_b = eps_b
        self.p = p
        self.d_l = d_l
        self.radiative = True
        self.M0 = E0 / (g0 * c**2)
        self.L = ((17 * self.M0) / (16 * pi * m_p * self.n)) ** (1 / 3)

    def gamma(self, t):
        # return 1000
        if self.radiative:
            return ((4 * c * t) / self.L) ** (-3 / 7)
        else:
            return (
                (17 * self.E(t)) / (1024 * pi * self.n * m_p * c**5 * t**3)
            ) ** (1 / 8)

    def B(self, t):
        return (32 * pi * m_p * self.eps_b * self.n) ** (1 / 2) * self.gamma(t) * c

    def g_m(self, t):
        return self.eps_e * ((self.p - 2) / (self.p - 1)) * (m_p / m_e) * self.gamma(t)

    def g_c(self, t):
        return (6 * pi * m_e * c) / (sigma_T * self.gamma(t) * self.B(t) ** 2 * t)

    def max_pow(self, t):
        return ((m_e * c**2 * sigma_T) / (3 * q_e)) * self.gamma(t) * self.B(t)

    def max_flux(self, t):
        return (self.N_e(t) * self.max_pow(t)) / (4 * pi * self.d_l**2)

    def frequency(self, t, gamma):
        return self.gamma(t) * gamma**2 * ((q_e * self.B(t)) / (2 * pi * m_e * c))

    def radius(self, t):
        if self.radiative:
            return ((4 * c * t) / self.L) ** (1 / 7) * self.L
        else:
            return ((17 * self.E(t) * t) / (4 * pi * m_p * self.n * c)) ** (1 / 4)

    def E(self, t):
        return (
            16 * pi * self.gamma(t) ** 2 * self.radius(t) ** 3 * self.n * m_p * c**2
        ) / 17

    def M(self, t):
        return self.E(t) / self.g0 * c**2

    def N_e(self, t):
        return (4 * pi * self.radius(t) ** 3 * self.n) / 3

    def calc_spectrum(self, t, freq):
        self.check_radiative(t)
        print(self.is_radiative)
        is_fast = self.g_m(t) > self.g_c(t)
        print(is_fast)
        spec = np.zeros_like(freq)

        spec[is_fast] = self.fast_cooling(t, freq[is_fast])
        spec[~is_fast] = self.slow_cooling(t, freq[~is_fast])

        return spec
        # if self.g_m(t) > self.g_c(t):
        #     return self.fast_cooling(t, freq)
        # else:
        #     return self.slow_cooling(t, freq)

    def slow_cooling(self, t, freq):
        crit = self.frequency(t, self.g_c(t))
        mini = self.frequency(t, self.g_m(t))

        flux = np.zeros_like(freq)

        min_gt_freq = mini > freq
        crit_gt_freq_a_freq_gt_mini = (crit > freq) & (freq > mini)
        freq_gt_crit = freq > crit

        flux[min_gt_freq] = (freq[min_gt_freq] / mini) ** (1 / 3) * self.max_flux(t)
        flux[crit_gt_freq_a_freq_gt_mini] = (
            freq[crit_gt_freq_a_freq_gt_mini] / mini
        ) ** (-(self.p - 1) / 2) * self.max_flux(t)
        flux[freq_gt_crit] = (
            (crit / mini) ** (-(self.p - 1) / 2)
            * (freq[freq_gt_crit] / crit) ** (-self.p / 2)
            * self.max_flux(t)
        )

        return flux
        # if mini > freq:
        #     return (freq / mini)**(1/3) * self.max_flux(t)
        # elif crit > freq and freq > mini:
        #     return (freq / mini)**(-(self.p-1)/2) * self.max_flux(t)
        # elif freq > crit:
        #     return (crit / mini)**(-(self.p-1)/2) * (freq / crit)**(-self.p / 2) * self.max_flux(t)

    def fast_cooling(self, t, freq):
        crit = self.frequency(t, self.g_c(t))
        mini = self.frequency(t, self.g_m(t))

        flux = np.zeros_like(freq)

        crit_gt_freq = crit > freq
        min_gt_freq_a_freq_gt_crit = (mini > freq) & (freq > crit)
        freq_gt_min = freq > mini

        flux[crit_gt_freq] = (freq[crit_gt_freq] / crit) ** (1 / 3) * self.max_flux(t)
        flux[min_gt_freq_a_freq_gt_crit] = (
            freq[min_gt_freq_a_freq_gt_crit] / crit
        ) ** (-1 / 2) * self.max_flux(t)
        flux[freq_gt_min] = (
            (mini / crit) ** (-1 / 2)
            * (freq[freq_gt_min] / mini) ** (-self.p / 2)
            * self.max_flux(t)
        )

        return flux
        # if crit > freq:
        #     return (freq / crit)**(1/3) * self.max_flux(t)
        # elif min > freq and freq > crit:
        #     return (freq / crit)**(-1/2) * self.max_flux(t)
        # elif freq > min:
        #     return (min / crit)**(-1/2) * (freq / min)**(-self.p / 2) * self.max_flux(t)

    def check_radiative(self, t):
        # self.is_radiative =  (self.g_c(t) < self.g_m(t)) and (self.eps_e >= 0.6)
        self.is_radiative = True