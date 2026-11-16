"""Figura 4 do TCC (Seção 2.3.2) – Esquema dos termos de 1ª e 2ª ordem da série de Born. (a) espalhamento único em r'; (b) espalhamento em r'', propagação livre G(r', r'')
e novo espalhamento em r'. Figura original do autor."""
import numpy as np
from estilo import *
from matplotlib.patches import Circle, FancyArrowPatch

fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
for ax in axs:
    ax.set_aspect("equal"); ax.axis("off"); ax.set(xlim=(-3.3, 3.4), ylim=(-2.1, 2.3))
    # região do potencial: disco sombreado com borda difusa
    for rr, al in zip(np.linspace(1.5, 0.2, 14), np.linspace(0.03, 0.10, 14)):
        ax.add_patch(Circle((0, 0), rr, color=AZUL, alpha=al, lw=0))
    ax.text(1.05, -1.45, r"$U(\mathbf{r})$", color=AZUL, fontsize=9)
    # onda plana incidente: frentes de onda
    for x in (-3.0, -2.7, -2.4, -2.1):
        ax.plot([x, x], [-1.2, 1.2], color=CINZA, lw=0.9)
    ax.annotate("", xy=(-1.65, 1.65), xytext=(-3.0, 1.65), arrowprops=dict(arrowstyle="-|>", color=TINTA, lw=1.0))
    ax.text(-2.35, 1.8, r"$\mathbf{k}_i$", ha="center", fontsize=10)
    ax.text(-2.55, -1.55, r"$\psi_{\mathrm{inc}}$", ha="center", fontsize=9, color=CINZA)
def ponto(ax, p, lab, dx=0.12, dy=0.12):
    ax.plot(*p, "o", color=LARANJA, ms=5.5, zorder=5); ax.text(p[0] + dx, p[1] + dy, lab, fontsize=10, color=LARANJA)
def seta(ax, a, b, ls="-", cor=TINTA, lw=1.1):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=10, color=cor, lw=lw, ls=ls, shrinkA=4, shrinkB=4))
# (a) primeira ordem
a1 = axs[0]; p1 = np.array([0.25, 0.35])
seta(a1, (-2.0, p1[1]), p1); ponto(a1, p1, r"$\mathbf{r}'$", dx=0.05, dy=-0.42)
d = np.array([np.cos(np.radians(38)), np.sin(np.radians(38))])
seta(a1, p1, p1 + 2.6 * d); a1.text(*(p1 + 2.65 * d + [0.05, 0.05]), r"$\mathbf{k}_f$", fontsize=10)
a1.text(0, -1.95, r"1ª ordem: $\int G\,U\,\psi_{\mathrm{inc}}$", ha="center", fontsize=9)
# (b) segunda ordem
a2 = axs[1]; p2 = np.array([-0.55, -0.45]); p1 = np.array([0.55, 0.55])
seta(a2, (-2.0, p2[1]), p2); ponto(a2, p2, r"$\mathbf{r}''$", dx=-0.15, dy=-0.38)
seta(a2, p2, p1, ls=(0, (3, 2)), cor=VERDE, lw=1.2); ponto(a2, p1, r"$\mathbf{r}'$", dx=-0.45, dy=0.12)
m = (p1 + p2) / 2; a2.text(m[0] + 0.12, m[1] - 0.2, r"$G(\mathbf{r}',\mathbf{r}'')$", fontsize=8.5, color=VERDE)
seta(a2, p1, p1 + 2.4 * d); a2.text(*(p1 + 2.45 * d + [0.05, 0.05]), r"$\mathbf{k}_f$", fontsize=10)
a2.text(0, -1.95, r"2ª ordem: $\iint G\,U\,G\,U\,\psi_{\mathrm{inc}}$", ha="center", fontsize=9)
rotulo(axs[0], "a"); rotulo(axs[1], "b")
fig.tight_layout()
salvar(fig, "2/fig_serie_born")
