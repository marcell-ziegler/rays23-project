from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt
from afterglow import Afterglow

m_e = const.m_e.to("g").value  # type:ignore
m_p = const.m_p.to("g").value  # type:ignore
pi = np.pi
c = const.c.to("cm/s").value  # type:ignore
sigma_T = const.sigma_T.cgs.value  # type: ignore
q_e = const.e.esu.value  # type: ignore


ag = Afterglow(E0=1e50, n=1e-3, g0=1000, eps_e=0.1, eps_b=0.2, p=2.5, d_l=1e28)

freqs = np.logspace(7, 30, base=10)
time = 1
densities = ag.calc_spectrum(time, freqs)

plt.axvline(ag.frequency(time, ag.g_m(time)), label="g_m")
plt.axvline(ag.frequency(time, ag.g_c(time)), ls="--", label="g_c")
nu_13 = freqs ** (1 / 3) / 1e32
nu_12 = freqs ** (-1 / 2) / 1e20
nu_p2 = freqs ** (-ag.p / 2) / 1
nu_mp12 = freqs ** (-(ag.p - 1) / 2) / 1e15
with open("data.csv", "w") as file:
    for line in freqs:
        file.write(str(line) + ",")
    file.write("\n\n\n")
    for line in densities:
        file.write(str(line) + ",")

print(ag.frequency(time, ag.g_c(time)), ag.frequency(time, ag.g_m(time)))
plt.loglog(freqs, densities)
plt.loglog(freqs, nu_13, label="nu^(1/3)")
plt.loglog(freqs, nu_12, label="nu^(-1/2)")
plt.loglog(freqs, nu_p2, label="nu^(-p/2)")
plt.loglog(freqs, nu_mp12, label="nu^(-(p-1)/2)")
# print(freqs, densities, sep="\n\n")
plt.legend()
plt.savefig("out.jpg")
