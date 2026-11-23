"""Figura 3 do TCC (Seção 2.3.1). O ouro está na energia das partículas de Bieler (7,35 MeV), que o usou como padrão – Regimes do espalhamento alfa-núcleo (parâmetro de Sommerfeld).

(a) diagrama D x lambda_bar_cm (log-log): eta = D/(2 lambda_bar_cm) -> isolinhas de inclinação +1;
    cada sistema percorre uma reta de inclinação +2 quando E varia (D ~ 1/E, lambda_bar ~ 1/sqrt(E)).
(b) D/R_toque x E_lab.

eta = k_e z Z e^2/(hbar v)          (v = velocidade relativa = velocidade da alfa no laboratório)
D   = k_e z Z e^2 / E_cm            (máxima aproximação, colisão frontal; E_cm ~ E para alvo pesado)
R_toque = 1,2 (A_alpha^{1/3} + A^{1/3}) fm  -> barreira coulombiana ~ energia em que D = R_toque
Fronteiras (convenção do texto): eta < 1 perturbativo coulombiano (Born para o Coulomb);
1 < eta < 4 semiclássico; eta > 4 clássico.
"""
import numpy as np
from estilo import *   # define o backend (pgf) antes de qualquer outro import do matplotlib
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from scatter_core import E2, ALPHA, M_ALPHA, HBARC

M_P = 938.272
alvos = [("Au", 79, 197, 196.967 * 931.494), ("Al", 13, 27, 26.982 * 931.494), ("H", 1, 1, M_P)]
cores = {"Au": AZUL, "Al": LARANJA, "H": VERDE}
estilos = {"Au": "-", "Al": "--", "H": "-."}
# experimentos de referência: (alvo, E_lab, rótulo curto, rótulo longo)
exps = [("Au", 7.35, "Bieler (Au)", "Bieler (1924), padrão\nAu, 7,35 MeV"),
        ("Al", 7.35, "Bieler (Al)", "Bieler (1924)\nAl, 7,35 MeV"),
        ("H", 7.68, "Chadwick–Bieler", "Chadwick e Bieler (1921)\nH, 7,68 MeV")]
E = np.logspace(np.log10(0.5), np.log10(50), 400)   # E_lab da alfa (MeV)
eta_de = lambda E, Z: 2 * Z * ALPHA / np.sqrt(2 * E / M_ALPHA)
dados = {nome: (Z, A, M) for nome, Z, A, M in alvos}
D_de = lambda E, nome: 2 * dados[nome][0] * E2 / (E * dados[nome][2] / (dados[nome][2] + M_ALPHA))
Rt_de = lambda nome: 1.2 * (4 ** (1 / 3) + dados[nome][1] ** (1 / 3))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.9))

# ---------------- (a) diagrama D x lambda_bar_cm: eta = D/(2 lambda_bar_cm) vira uma família de retas de inclinação +1
mu_de = lambda nome: M_ALPHA * dados[nome][2] / (M_ALPHA + dados[nome][2])
lb_de = lambda E, nome: HBARC / np.sqrt(2 * mu_de(nome) * E * dados[nome][2] / (dados[nome][2] + M_ALPHA))
xmin, xmax, ymin, ymax = 0.2, 30, 0.2, 800
xx = np.logspace(np.log10(xmin), np.log10(xmax), 200)
a1.fill_between(xx, 2 * xx, ymax, color="#eda100", alpha=0.10, lw=0)   # eta > 1
a1.fill_between(xx, ymin, 2 * xx, color=AZUL, alpha=0.09, lw=0)       # eta < 1
for et in (0.1, 10, 100):
    a1.plot(xx, 2 * et * xx, color=CINZA, lw=0.5, ls=(0, (1, 2)))
a1.plot(xx, 8 * xx, color=CINZA, lw=0.8, ls=":")
a1.plot(xx, 2 * xx, color=TINTA, lw=1.3, ls=(0, (5, 2.5)), zorder=3)
a1.set(xscale="log", yscale="log", xlim=(xmin, xmax), ylim=(ymin, ymax),
       xlabel=r"$\bar{\lambda}_\mathrm{cm} = \hbar/\sqrt{2\mu E_\mathrm{cm}}$ [fm]", ylabel=r"$D = k_e\, zZe^2/E_\mathrm{cm}$ [fm]")
fig.canvas.draw()
def rotulo_isolinha(ax, eta, x, texto, **kw):
    """Escreve o texto sobre a reta D = 2 eta x, girado com a inclinação aparente na tela."""
    p1 = ax.transData.transform((x, 2 * eta * x)); p2 = ax.transData.transform((x * 1.5, 2 * eta * x * 1.5))
    ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    ax.text(x, 2 * eta * x * 1.18, texto, rotation=ang, rotation_mode="anchor", ha="left", va="bottom", **kw)
rotulo_isolinha(a1, 1, 6.0, r"$\eta = 1$", fontsize=8, color=TINTA)
rotulo_isolinha(a1, 4, 6.0, r"$\eta = 4$", fontsize=7.5, color=CINZA)
for et, x in [(0.1, 12.0), (10, 0.25), (100, 0.25)]:
    rotulo_isolinha(a1, et, x, (f"$\\eta = {et:g}$").replace(".", "{,}"), fontsize=6.5, color=CINZA)
a1.text(0.03, 0.97, "clássico / semiclássico\n" + r"($\eta > 1$)", transform=a1.transAxes, fontsize=7.5, va="top")
a1.text(0.97, 0.03, "perturbativo coulombiano\n" + r"Born ($\eta < 1$)", transform=a1.transAxes, fontsize=7.5, va="bottom", ha="right")
for nome, Z, A, M in alvos:
    lb = lb_de(E, nome); D = D_de(E, nome)
    a1.plot(lb, D, color=cores[nome], ls=estilos[nome])
    inc = np.polyfit(np.log(lb), np.log(D), 1)[0]
    print(f"{nome}: inclinação da curva ln D x ln lambda_bar = {inc:.4f} (esperado 2)")
pos_txt = {"Au": ((1.15, 50), "left", "center"), "Al": ((0.85, 6.6), "right", "bottom"), "H": ((4.7, 1.75), "left", "top")}
for nome, x, curto, _ in exps:
    lb = lb_de(x, nome); D = D_de(x, nome)
    a1.plot([lb], [D], "o", ms=6, color=TINTA, mec="white", mew=1, zorder=5)
    (tx_, ty_), ha_, va_ = pos_txt[nome]
    a1.text(tx_, ty_, curto, fontsize=7, ha=ha_, va=va_, zorder=6,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
    print(f"{curto}: lambda_bar_cm = {lb:.3f} fm, D = {D:.3f} fm, eta = D/(2 lambda_bar) = {D/(2*lb):.3f}"
          f"  (eta direto = {eta_de(x, dados[nome][0]):.3f})")
# Seção 4.1.1: ouro e alumínio na mesma energia (a das partículas de Bieler no alumínio)
print(f"7,35 MeV: beta = {np.sqrt(2 * 7.35 / M_ALPHA):.4f}, eta(Au) = {eta_de(7.35, 79):.2f}, eta(Al) = {eta_de(7.35, 13):.2f}")
# checagem: a reta eta = 1 passa por (1 fm, 2 fm) e tem inclinação 1
print("isolinha eta=1 em lambda=1 fm:", 2 * 1 * 1.0, "fm; inclinação:", np.polyfit(np.log(xx), np.log(2 * xx), 1)[0])
# seta: sentido de energia crescente
a1.annotate("", xy=(0.4, 17), xytext=(0.8, 68), arrowprops=dict(arrowstyle="->", color=CINZA, lw=0.8))
a1.text(0.8, 80, r"$E$ cresce", fontsize=7, color=CINZA, ha="center", va="bottom")

# ---------------- (b) D / R_toque x E: uma única linha de referência (D = R_toque)
ymin_b, ymax_b = 0.05, 50
a2.axhspan(ymin_b, 1, facecolor="none", edgecolor=CINZA, hatch="////", lw=0, alpha=0.35)
a2.axhline(1, color=TINTA, lw=1.3, ls=(0, (5, 2.5)))
a2.text(0.03, 0.12, "$D \\leq R_\\mathrm{toque}$: contato nuclear\nem colisão frontal", transform=a2.get_yaxis_transform(),
        fontsize=7.5, va="center", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
a2.text(0.03, 1.12, "$D = R_\\mathrm{toque}$", transform=a2.get_yaxis_transform(), fontsize=8, va="bottom")
linhas = []
for nome, Z, A, M in alvos:
    Rt = Rt_de(nome)
    exp = [e for e in exps if e[0] == nome][0]
    Elab = exp[1]; Ecm = Elab * M / (M + M_ALPHA)
    lab = (fr"$\alpha + {{}}^{{{A}}}\mathrm{{{nome}}}$ ($R_\mathrm{{toque}} = " + f"{Rt:.1f}".replace(".", "{,}") + r"$ fm) — "
           + exp[2] + ": " + (f"$E_\\mathrm{{lab}} = {Elab:.2f}$ MeV, $E_\\mathrm{{cm}} = {Ecm:.2f}$ MeV").replace(".", "{,}"))
    l, = a2.plot(E, D_de(E, nome) / Rt, color=cores[nome], ls=estilos[nome], label=lab)
    linhas.append(l)
offs_b = {"Au": (6, 4), "Al": (6, 4), "H": (-8, -2)}
for nome, x, curto, _ in exps:
    y = D_de(x, nome) / Rt_de(nome)
    a2.plot([x], [y], "o", ms=6, color=TINTA, mec="white", mew=1, zorder=5)
    o = offs_b[nome]
    M = dados[nome][2]; Ecm = x * M / (M + M_ALPHA)
    a2.annotate(curto, (x, y), xytext=o, textcoords="offset points", fontsize=7,
                ha="right" if o[0] < 0 else "left", va="center" if o[0] < 0 else "bottom",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=6)
    print(f"{curto}: E_cm = {Ecm:.3f} MeV, D = {D_de(x, nome):.2f} fm, R_toque = {Rt_de(nome):.2f} fm, D/R_toque = {y:.2f}")
a2.set(xscale="log", yscale="log", ylim=(ymin_b, ymax_b),
       xlabel=r"$E_\alpha$ (lab) [MeV]", ylabel=r"$D/R_\mathrm{toque}$")

a2.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50]))
a2.xaxis.set_major_formatter(FixedFormatter(["0,5", "1", "2", "5", "10", "20", "50"]))
a1.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5, 10, 20]))
a1.xaxis.set_major_formatter(FixedFormatter(["0,2", "0,5", "1", "2", "5", "10", "20"]))
a1.yaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]))
a1.yaxis.set_major_formatter(FixedFormatter(["0,2", "0,5", "1", "2", "5", "10", "20", "50", "100", "200", "500"]))
for ax in (a1, a2):
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
a2.yaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]))
a2.yaxis.set_major_formatter(FixedFormatter(["0,05", "0,1", "0,2", "0,5", "1", "2", "5", "10", "20", "50"]))
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout(rect=(0, 0.13, 1, 1))
fig.legend(handles=linhas, loc="lower center", ncol=1, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
salvar(fig, "4.1/fig01_regimes_sommerfeld_v2")
