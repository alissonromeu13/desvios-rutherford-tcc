"""Figura 9 do TCC (Seção 2.4.5) – Fatores de forma da Seção 2.4 comparados a uma densidade realista (Fermi de 2 parâmetros) com o
MESMO raio quadrático médio. Ilustração para o chumbo: c = 6,62 fm, z = 0,546 fm, parâmetros 2pF do 207Pb
(De Vries, De Jager & De Vries, 1987, Tab. I, p. 510; o 208Pb só tem análises FB/SOG). F(q) = (4 pi/q) int rho(r) r sin(qr) dr, normalizado.
"""
import numpy as np
from scipy.integrate import simpson
from estilo import *

r = np.linspace(1e-6, 25, 20001)
c, z = 6.62, 0.546   # 2pF do 207Pb (De Vries et al. 1987, Tab. I, p. 510); o 208Pb só tem análises FB/SOG
rho_f = 1 / (1 + np.exp((r - c) / z))
rms = np.sqrt(simpson(rho_f * r**4, x=r) / simpson(rho_f * r**2, x=r))
R_hs = np.sqrt(5 / 3) * rms           # esfera homogênea com mesmo <r^2>
a_g = rms * np.sqrt(2 / 3)            # gaussiana exp(-r^2/a^2): <r^2> = 3a^2/2
a_e = rms / np.sqrt(12)               # exponencial exp(-r/a): <r^2> = 12 a^2
print(f"r_rms(Fermi) = {rms:.3f} fm -> R_esfera = {R_hs:.3f} fm, a_gauss = {a_g:.3f} fm, a_exp = {a_e:.3f} fm")

q = np.linspace(1e-3, 3.0, 1500)
def F_num(rho):
    norm = simpson(rho * r**2, x=r)
    return np.array([simpson(rho * r * np.sin(qq * r), x=r) / qq for qq in q]) / norm
F_f = F_num(rho_f)
x = q * R_hs; F_hs = 3 * (np.sin(x) - x * np.cos(x)) / x**3
F_g = np.exp(-(q * a_g) ** 2 / 4)
F_e = 1 / (1 + (q * a_e) ** 2) ** 2
# checagem: F numérico da esfera = analítico
rho_hs = (r < R_hs).astype(float)
print("checagem esfera (numérico vs analítico) max dif =", np.max(np.abs(F_num(rho_hs) - F_hs)))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw={"width_ratios": [1, 1.5]})
a1.plot(r, rho_f, color=TINTA, label="Fermi (2pF)")
a1.plot(r, rho_hs * simpson(rho_f * r**2, x=r) / simpson(rho_hs * r**2, x=r), color=AZUL, ls="--", label="esfera homogênea (Seç. 2.4.1)")
N_f = simpson(rho_f * r**2, x=r)
re_ = np.exp(-r / a_e); re_n = re_ * N_f / simpson(re_ * r**2, x=r)
rg = np.exp(-(r / a_g) ** 2); rg_n = rg * N_f / simpson(rg * r**2, x=r)
print(f"rho(0)/rho0 com o mesmo <r^2> e a mesma carga: exponencial = {re_n[0]:.2f}, gaussiana = {rg_n[0]:.2f}, esfera = {N_f/(R_hs**3/3):.3f}")
a1.plot(r, re_n, color=LARANJA, ls=":", label="exponencial/dipolo (Seç. 2.4.2)")
a1.plot(r, rg_n, color=VERDE, ls="-.", label="gaussiana (Seç. 2.4.3)")
# escala log: a exponencial, com o mesmo <r^2>, tem rho(0) ~ 13 rho0 e não cabe na escala linear
a1.set(yscale="log", xlabel=r"$r$ [fm]", ylabel=r"$\rho(r)/\rho_0$", xlim=(0, 12), ylim=(1e-3, 30))
for F, cc, ls in [(F_f, TINTA, "-"), (F_hs, AZUL, "--"), (F_e, LARANJA, ":"), (F_g, VERDE, "-.")]:
    a2.plot(q, F**2, color=cc, ls=ls)
a2.set(yscale="log", ylim=(1e-8, 1.5), xlim=(0, 3), xlabel=r"$q$ [fm$^{-1}$]", ylabel=r"$|F(q)|^2$")
virgula(a1, "x"); virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
h1, l1 = a1.get_legend_handles_labels()
fig.tight_layout(rect=(0, 0.11, 1, 1))
fig.legend(h1, l1, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=7.5, frameon=False,
           title=r"$^{208}$Pb (2pF do $^{207}$Pb), mesma carga e mesmo $\langle r^2\rangle^{1/2}$ = " + f"{rms:.2f}".replace(".", "{,}") + " fm",
           title_fontsize=7.5, handlelength=2.6)
salvar(fig, "5.4/fig10_fatores_forma_realistas_v4")
