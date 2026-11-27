"""Figura 12 do TCC (Seção 4.2.4.2) – Disco negro: S_l = 0 para l + 1/2 < kR, S_l = 1 caso contrário.
 Ondas parciais exatas: f = (i/2k) sum (2l+1) P_l(cos th),  sigma_el = sigma_r = pi L^2/k^2.
 Eiconal: f = i k R^2 J1(qR)/(qR), q = k theta (P_l(cos th) ~ J0((l+1/2) th) para pequenos ângulos, Seção 4.2.4)
   -> (dsigma/dOmega)/R^2 = (kR)^2 [J1(k theta R)/(k theta R)]^2  (padrão de Airy).
   Com q = k theta os mínimos coincidem com os da soma de ondas parciais; com q = 2k sen(theta/2) eles se deslocam.
 Painel único (v1.2): o painel das seções de choque integradas, com borda suavizada, foi retirado por
 exceder o que a seção deduz; sigma_el = sigma_r = pi R^2 é obtido analiticamente no texto.
 Teorema óptico: sigma_tot = (4pi/k) Im f(0) = 2 pi R^2 para kR >> 1.
"""
import numpy as np
from scipy.special import j1, eval_legendre
from estilo import *

fig, a1 = plt.subplots(figsize=(4.6, 3.2))
th = np.radians(np.linspace(0.05, 45, 900))
for kR, c in [(10, AZUL), (30, LARANJA)]:
    k, R = kR, 1.0
    L = int(np.floor(kR + 0.5))
    f = sum((2 * l + 1) * eval_legendre(l, np.cos(th)) for l in range(L)) * 1j / (2 * k)
    q = k * th                                            # aproximação de pequenos ângulos
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
a1.set(yscale="log", ylim=(1e-3, 3e2), xlim=(0, 45), xlabel=r"$\theta$ [graus]", ylabel=r"$(\mathrm{d}\sigma/\mathrm{d}\Omega)/R^2$")
a1.legend(fontsize=7, loc="upper right")

virgula(a1, "x")
fig.tight_layout()
salvar(fig, "6.3/fig05_disco_negro_difracao_v3")
