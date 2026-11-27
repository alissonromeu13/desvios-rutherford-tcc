"""Figura 7 do TCC (Seção 2.4.3) – Distribuição gaussiana rho ~ exp(-r²/a²) e fator de forma F(x) = exp(-x²/4), x = qa.
<r²> = 3a²/2  =>  F ≈ 1 − q²<r²>/6 = 1 − x²/4 para q pequeno.
(a) F e |F|² em escala linear, com a aproximação de q pequeno; (b) ln|F|² = −x²/2: reta em função de x²
(decaimento mais rápido que qualquer potência de q). Cor: verde = gaussiana (como na Figura 8)."""
import numpy as np
from estilo import *

x = np.linspace(0, 6, 2000)
F = np.exp(-x**2 / 4)
print("F = 1/2 em x =", round(2 * np.sqrt(np.log(2)), 3), "; |F|^2 = 1e-4 em x =", round(np.sqrt(2 * np.log(1e4)), 2))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
a1.plot(x, F, color=VERDE, label=r"$F(x)$")
a1.plot(x, F**2, color=VERDE, ls="-.", lw=1.5, label=r"$|F(x)|^2$")
xs = np.linspace(0, 1.6, 100)
h_q, = a1.plot(xs, 1 - xs**2 / 4, color=CINZA, lw=1.0, ls=(0, (5, 2)), label=r"$q$ pequeno: $1 - q^2\langle r^2\rangle/6 = 1 - x^2/4$")
a1.set(xlim=(0, 5), ylim=(0, 1.05), xlabel=r"$x = qa$", ylabel="fator de forma")
x2 = np.linspace(0, 6, 400)
a2.plot(x2**2, np.exp(-x2**2 / 2), color=VERDE, ls="-.", lw=1.5)
a2.set(yscale="log", xlim=(0, 36), ylim=(1e-8, 1.5), xlabel=r"$x^2 = (qa)^2$", ylabel=r"$|F(x)|^2$")
a2.set_xticks([0, 6, 12, 18, 24, 30, 36])
a2.text(0.95, 0.92, r"$|F|^2 = e^{-x^2/2}$: reta em $x^2$", transform=a2.transAxes, ha="right", va="top", fontsize=7.5)
sec = a2.secondary_xaxis("top", functions=(np.sqrt, np.square))
sec.set_xticks([2, 3, 4, 5, 6]); sec.set_xlabel(r"$x = qa$", fontsize=8); sec.tick_params(labelsize=7.5)
fig.tight_layout(rect=(0, 0.09, 1, 1))
h = a1.get_legend_handles_labels()[0][:2] + [h_q]
fig.legend(handles=h, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, 0.0), handlelength=2.4)
virgula(a1); virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
salvar(fig, "5.3/fig_5_3_gaussiana")
