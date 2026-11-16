"""Figura 5 do TCC (Seção 2.4.1) – Fator de forma da esfera homogênea, F(x) = 3 j1(x)/x = 3 (sin x − x cos x)/x³, x = qR.
(a) F e |F|² em escala linear; (b) |F|² em escala log, com a envoltória assintótica 9/x⁴
(j1(x) ≈ −cos x / x para x ≫ 1  =>  F ≈ −3 cos x / x²). Zeros: tan x = x."""
import numpy as np
from scipy.optimize import brentq
from estilo import *

x = np.linspace(1e-4, 15, 3000)
F = 3 * (np.sin(x) - x * np.cos(x)) / x**3
zeros = [brentq(lambda t: np.tan(t) - t, (n + 0.5) * np.pi - 1.2, (n + 0.5) * np.pi - 1e-6) for n in range(1, 5)]
print("zeros de j1:", np.round(zeros, 3))
xm = [brentq(lambda t: 3*((t**2-3)*np.sin(t)+3*t*np.cos(t))/t**4, a, b) for a, b in [(5, 6.5), (8, 9.8), (11, 13)]]
Fm = [3*(np.sin(t)-t*np.cos(t))/t**3 for t in xm]
print("extremos de F:", np.round(xm, 3), np.round(Fm, 4), "|F|^2:", np.round(np.square(Fm), 5))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
a1.plot(x, F, color=AZUL, label=r"$F(x)$")
a1.plot(x, F**2, color=AZUL, ls="--", lw=1.4, label=r"$|F(x)|^2$")
a1.axhline(0, color=CINZA, lw=0.7)
a2.plot(x, F**2, color=AZUL, ls="--", lw=1.4, label=r"$|F(x)|^2$")
x2 = np.linspace(3, 15, 400)
a2.plot(x2, 9 / x2**4, color=CINZA, ls=":", lw=1.2, label=r"envoltória $9/x^4$")
for ax in (a1, a2):
    for z0 in zeros:
        ax.axvline(z0, color=CINZA, lw=0.6, ls=(0, (1, 2)))
    ax.set(xlim=(0, 15), xlabel=r"$x = qR$")
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12, 14])
for z0 in zeros:
    a1.text(z0, 1.03, f"{z0:.2f}".replace(".", ","), fontsize=6.5, color=CINZA, ha="center", va="bottom")
a1.set(ylim=(-0.2, 1.1), ylabel="fator de forma")
a2.set(yscale="log", ylim=(1e-4, 1.5), ylabel=r"$|F(x)|^2$")
fig.tight_layout(rect=(0, 0.07, 1, 1))
h = a1.get_legend_handles_labels()[0][:1] + a2.get_legend_handles_labels()[0]
fig.legend(handles=h, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
virgula(a1); virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
salvar(fig, "5.1/fig_5_1_esfera_homogenea")
