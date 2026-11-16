"""Figura 13 do TCC (Seção 4.2.4.2) – Disco negro: S_l = 0 para l + 1/2 < kR, S_l = 1 caso contrário.
 Ondas parciais exatas: f = (i/2k) sum (2l+1) P_l(cos th),  sigma_el = sigma_r = pi L^2/k^2.
 Eiconal: f = i k R^2 J1(qR)/(qR), q = k theta (aproximação de Hilb P_l(cos th) ~ J0((l+1/2) th), como na Seção 4.2.4.2)
   -> (dsigma/dOmega)/R^2 = (kR)^2 [J1(k theta R)/(k theta R)]^2  (padrão de Airy).
   Com q = k theta os mínimos coincidem com os da soma de ondas parciais; com q = 2k sen(theta/2) eles se deslocam.
 Painel (b): o disco negro só faz sentido para lambda_bar << R; as curvas começam em kR = 2.
 Teorema óptico: sigma_tot = (4pi/k) Im f(0) = 2 pi R^2 para kR >> 1.
"""
import numpy as np
from scipy.special import j1, eval_legendre
from estilo import *

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
th = np.radians(np.linspace(0.05, 45, 900))
for kR, c in [(10, AZUL), (30, LARANJA)]:
    k, R = kR, 1.0
    L = int(np.floor(kR + 0.5))
    f = sum((2 * l + 1) * eval_legendre(l, np.cos(th)) for l in range(L)) * 1j / (2 * k)
    q = k * th                                            # aproximação de pequenos ângulos (Hilb)
    pw = np.abs(f) ** 2 / R**2
    eik = (k * R * j1(q * R) / (q * R)) ** 2              # (dsigma/dOmega)/R^2
    a1.plot(np.degrees(th), pw, color=c, lw=1.6, label=f"$kR={kR}$ – ondas parciais")
    a1.plot(np.degrees(th), eik, color=TINTA, lw=0.8, ls="--")
    m = np.degrees(th) < 30
    q2 = 2 * k * np.sin(th / 2); eik2 = (k * R * j1(q2 * R) / (q2 * R)) ** 2
    print(f"kR={kR}: erro rms em log(dsigma/dOmega), theta < 30 graus: q = k theta -> "
          f"{np.sqrt(np.mean(np.log(eik[m] / pw[m]) ** 2)):.3f};  q = 2k sen(theta/2) -> "
          f"{np.sqrt(np.mean(np.log(eik2[m] / pw[m]) ** 2)):.3f}")
    # optical theorem check
    f0 = 1j / (2 * k) * L**2
    print(f"kR={kR}: sigma_tot (teor. óptico)/piR^2 = {4*np.pi/k*f0.imag/np.pi:.3f};  2(L/kR)^2 = {2*(L/kR)**2:.3f}")
# primeiro zero difrativo: k theta_1 R = x_1 (1o zero de J_1)
from scipy.special import jn_zeros
x1 = jn_zeros(1, 1)[0]
for kR, c in [(10, AZUL), (30, LARANJA)]:
    t1 = np.degrees(x1 / kR)
    a1.plot([t1], [1.0e-3], marker="^", ms=6, color=c, mec="white", mew=0.6, clip_on=False, zorder=7)
    a1.text(t1 + 0.6, 1.6e-3, r"$\theta_1$", fontsize=7.5, color=c, va="bottom",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
    print(f"kR={kR}: theta_1 = {t1:.2f} graus (x_1 = {x1:.4f})")
a1.plot([], [], color=TINTA, lw=0.8, ls="--", label=r"eiconal: $(kR)^2[J_1(k\theta R)/k\theta R]^2$")
a1.set(yscale="log", ylim=(1e-3, 3e2), xlim=(0, 45), xlabel=r"$\theta$ [graus]", ylabel=r"$(d\sigma/d\Omega)/R^2$")
a1.legend(fontsize=7, loc="upper right")

KR_MIN = 2.0
kRs = np.linspace(KR_MIN, 40, 400)
Delta = 0.5   # borda suavizada: S_l = 1/(1+exp((kR - l - 1/2)/Delta)) remove o "dente de serra" do l discreto
se, sr = [], []
for x in kRs:
    l = np.arange(int(x + 15 * Delta + 5))
    S = 1 / (1 + np.exp((x - l - 0.5) / Delta))
    se.append(np.sum((2 * l + 1) * (1 - S) ** 2) / x**2)
    sr.append(np.sum((2 * l + 1) * (1 - S**2)) / x**2)
se, sr = np.array(se), np.array(sr)
a2.plot(kRs, se + sr, color=TINTA, lw=1.4, label=r"$\sigma_\mathrm{tot}$")
a2.plot(kRs, sr, color=AZUL, lw=1.4, ls="--", label=r"$\sigma_\mathrm{reação}$")
a2.plot(kRs, se, color=VERDE, lw=1.4, ls="-.", label=r"$\sigma_\mathrm{el}$ (sombra)")
print("kR=40: sigma_el, sigma_r, tot /piR^2 =", se[-1].round(3), sr[-1].round(3), (se+sr)[-1].round(3))
a2.axhline(2, color=LARANJA, lw=0.8, ls=":"); a2.axhline(1, color=CINZA, lw=0.8, ls=":")
a2.text(10, 2.08, r"$2\pi R^2$", color=LARANJA, fontsize=8, va="bottom")
a2.text(24, 0.75, r"$\pi R^2$ (área geométrica)", color=CINZA, fontsize=7.5)
a2.axvspan(0, KR_MIN, facecolor="none", edgecolor=CINZA, hatch="////", lw=0, alpha=0.4)
a2.text(KR_MIN + 0.6, 0.25, r"à esquerda: $\bar{\lambda} \gtrsim R/2$," + "\nmodelo não se aplica", fontsize=6.8, color=CINZA)
a2.set(xlabel=r"$kR$", ylabel=r"$\sigma/\pi R^2$", ylim=(0, 3.0), xlim=(0, 40))
a2.legend(fontsize=7.5, loc="upper right")
virgula(a1, "x"); virgula(a2)
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout()
salvar(fig, "6.3/fig05_disco_negro_difracao_v3")
