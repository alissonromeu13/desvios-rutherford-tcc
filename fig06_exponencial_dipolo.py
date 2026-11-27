"""Figura 6 do TCC (Seção 2.4.2) – Distribuição exponencial rho ~ exp(-r/a) e fator de forma de dipolo F(x) = 1/(1+x²)², x = qa.
(a) F e |F|² em escala linear, com a aproximação de q pequeno F ≈ 1 − q²<r²>/6 = 1 − 2x² (<r²> = 12 a²);
(b) |F|² em escala log, com a assíntota x^-8. Sem zeros. Cores: laranja = exponencial (como na Figura 8).
A comparação com os outros modelos, com o mesmo <r²>, fica na Seção 2.4.4 (Figura 8)."""
import numpy as np
from estilo import *

x = np.linspace(1e-4, 6, 2000)
F = 1 / (1 + x**2) ** 2
print("F = 1/2 em x =", round(np.sqrt(np.sqrt(2) - 1), 3), "; erro de 1-2x² em x=0,2:",
      round(abs((1 - 2 * 0.04) - 1 / 1.04**2) / (1 / 1.04**2), 4), "; |F|²/x^-8 em x=6:", round((1/(1+36)**4) * 6**8, 3))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
a1.plot(x, F, color=LARANJA, label=r"$F(x)$")
a1.plot(x, F**2, color=LARANJA, ls=":", lw=1.9, label=r"$|F(x)|^2$")
xs = np.linspace(0, 0.62, 100)
h_q, = a1.plot(xs, 1 - 2 * xs**2, color=CINZA, lw=1.0, ls=(0, (5, 2)), label=r"$q$ pequeno: $1 - q^2\langle r^2\rangle/6 = 1 - 2x^2$")
a1.set(xlim=(0, 4), ylim=(0, 1.05), xlabel=r"$x = qa$", ylabel="fator de forma")
a2.plot(x, F**2, color=LARANJA, ls=":", lw=1.9)
x2 = np.linspace(1.5, 6, 200)
h_a, = a2.plot(x2, x2**-8.0, color=CINZA, lw=1.0, ls=(0, (1, 1.5)), label=r"assíntota $x^{-8}$")
a2.set(yscale="log", xlim=(0, 6), ylim=(1e-7, 1.5), xlabel=r"$x = qa$", ylabel=r"$|F(x)|^2$")
fig.tight_layout(rect=(0, 0.09, 1, 1))
h = a1.get_legend_handles_labels()[0][:2] + [h_q, h_a]
fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, 0.0), handlelength=2.4)
virgula(a1); virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
salvar(fig, "5.2/fig_5_2_exponencial_dipolo")
