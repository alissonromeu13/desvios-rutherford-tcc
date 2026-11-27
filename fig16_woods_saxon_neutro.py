"""Figura 13 do TCC (Seção 4.4.2) – Partícula neutra no potencial óptico de Woods-Saxon: n + 208Pb, E_lab = 14 MeV.

Potencial óptico local de Koning & Delaroche, Nucl. Phys. A 713 (2003) 231, SEM o termo spin-órbita:
  U(r) = -V_V f(r;R_V,a_V) - i W_V f(r;R_V,a_V) + 4 i a_D W_D d/dr f(r;R_D,a_D)        (Eqs. 1-3, p. 240-241)
Profundidades pela Eq. 7 (p. 242), com E = energia de laboratório e E_f = E_f^n:
  V_V = v1[1 - v2 x + v3 x^2 - v4 x^3],  W_V = w1 x^2/(x^2 + w2^2),  W_D = d1 x^2/(x^2 + d3^2) exp(-d2 x),  x = E - E_f
Parâmetros do 208Pb: Tabela 5 (p. 246) e Tabela 6 (p. 247).
Numerov + casamento com Riccati-Hankel em dois pontos; checagens: teorema óptico e convergência.
"""
import numpy as np
from estilo import *
from scatter_core import smatrix, nuclear_amplitude, woods_saxon, HBARC, M_N, AMU

# ---------------------------------------------------------------- KD03, n + 208Pb
rV, aV, v1, v2, v3, v4, w1, w2, rD, aD = 1.244, 0.646, 50.6, 0.0069, 0.000015, 7e-9, 15.6, 88.0, 1.246, 0.510  # Tab. 5
d1, d2, d3, Efn = 13.8, 0.0180, 13.8, -5.65                                                                      # Tab. 6
A, E = 208, 14.0
x = E - Efn
VV = v1 * (1 - v2 * x + v3 * x**2 - v4 * x**3)
WV = w1 * x**2 / (x**2 + w2**2)
WD = d1 * x**2 / (x**2 + d3**2) * np.exp(-d2 * x)
RV, RD = rV * A ** (1 / 3), rD * A ** (1 / 3)
print(f"KD03 n+208Pb, E = {E} MeV (E - Ef = {x:.2f}): V_V = {VV:.2f}, W_V = {WV:.3f}, W_D = {WD:.2f} MeV; "
      f"R_V = {RV:.3f} fm, R_D = {RD:.3f} fm")

mu = M_N * A * AMU / (M_N + A * AMU)
Ecm = E * A * AMU / (A * AMU + M_N)

def U(r, absorve):
    xd = (r - RD) / aD
    g = np.exp(xd) / (1 + np.exp(xd)) ** 2            # g = -a_D df/dr
    W = (WV * woods_saxon(r, RV, aV) + 4 * WD * g) if absorve else 0.0
    return -VV * woods_saxon(r, RV, aV) - 1j * W

Rcorte = RV + 15 * aV
th = np.radians(np.linspace(0.5, 180, 700))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.5, 1]})
res = {}
for abs_, c, ls, nome in [(False, AZUL, "--", "real ($W = 0$)"), (True, LARANJA, "-", "complexo ($W_V$ e $W_D$)")]:
    lmax = 30
    S, k = smatrix(lambda r: U(r, abs_), mu, Ecm, lmax, Rcorte, h=0.01)
    S2, _ = smatrix(lambda r: U(r, abs_), mu, Ecm, lmax + 15, RV + 20 * aV, h=0.005)
    l = np.arange(lmax + 1); l2 = np.arange(len(S2))
    sel = np.pi / k**2 * np.sum((2 * l + 1) * np.abs(1 - S) ** 2) * 10
    sr = np.pi / k**2 * np.sum((2 * l + 1) * (1 - np.abs(S) ** 2)) * 10
    sel2 = np.pi / k**2 * np.sum((2 * l2 + 1) * np.abs(1 - S2) ** 2) * 10
    sr2 = np.pi / k**2 * np.sum((2 * l2 + 1) * (1 - np.abs(S2) ** 2)) * 10
    f0 = nuclear_amplitude(S, 0.0, k, np.array([0.0]))[0]
    sot = 4 * np.pi / k * f0.imag * 10
    conv = max(abs(sel2 - sel) / sel, abs(sr2 - sr) / max(sr, 1e-30) if abs_ else 0.0)
    assert abs(sot - (sel + sr)) / sot < 1e-8 and conv < 1e-4
    print(f"{nome}: sigma_el={sel:.0f} mb  sigma_r={sr:.0f} mb  tot={sel+sr:.0f} mb  teor.óptico={sot:.0f} mb  "
          f"(k={k:.4f}/fm, kR_V={k*RV:.2f}, kR_corte={k*Rcorte:.1f}); convergência relativa: {conv:.1e}")
    res[abs_] = (S, sel, sr)
    f = nuclear_amplitude(S, 0.0, k, th)
    a1.plot(np.degrees(th), np.abs(f) ** 2 * 10, color=c, ls=ls, label=nome)
a1.set(yscale="log", xlabel=r"$\theta_\mathrm{cm}$ [graus]", ylabel=r"$\mathrm{d}\sigma/\mathrm{d}\Omega$ [mb/sr]", xlim=(0, 180), ylim=(1e-2, 4e4))
a1.set_xticks([0, 30, 60, 90, 120, 150, 180])
a1.legend(fontsize=7.5, loc="upper right", title=r"n + $^{208}$Pb, $E = 14$ MeV", title_fontsize=7.5)
fmt = lambda v: f"{v/1000:.2f}".replace(".", "{,}")
_, se0, _ = res[False]; _, se1, sr1 = res[True]
tab = (r"\begin{tabular}{@{}lccc@{}} & $\sigma^{\mathrm{elás}}$ & $\sigma^{\mathrm{reação}}$ & $\sigma^{\mathrm{tot}}$ \\[1pt]"
       + r"real & $" + fmt(se0) + r"$ & 0 & $" + fmt(se0) + r"$ \\"
       + r"complexo & $" + fmt(se1) + r"$ & $" + fmt(sr1) + r"$ & $" + fmt(se1 + sr1) + r"$ \\ \end{tabular}")
a1.text(0.84, 0.03, tab, transform=a1.transAxes, fontsize=6.8, va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=CINZA, lw=0.5))
a1.text(0.84, 0.205, "seções de choque integradas [b]:", transform=a1.transAxes, fontsize=6.8, va="bottom", ha="right")
for abs_, c, ls in [(False, AZUL, "--"), (True, LARANJA, "-")]:
    S = res[abs_][0]
    a2.plot(np.arange(len(S)), np.abs(S), color=c, ls=ls, marker="o", ms=3)
eta = np.abs(res[True][0]); lmin = int(np.argmin(eta[:15]))
a2.set(xlabel=r"$l$", ylabel=r"$|S_l|$", ylim=(0, 1.05), xlim=(0, 20))
a2.axvline(k * RV, color=TINTA, lw=1.0, ls=(0, (4, 2)))
a2.text(k * RV + 0.35, 0.03, r"$l \approx kR \approx " + f"{k*RV:.1f}".replace(".", "{,}") + "$", fontsize=7.5)
a2.annotate("absorção máxima\n" + r"em $l \approx kR$", xy=(lmin, eta[lmin]), xytext=(10.5, 0.30),
            fontsize=6.8, va="center", arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA))
a2.text(19.6, 0.66, r"$|S_l| \to 1$ para ondas" + "\nque não alcançam o núcleo", fontsize=6.8, ha="right", va="top", color=CINZA)
virgula(a2, "y")
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout()
salvar(fig, "6.4/fig07_woods_saxon_neutro_v5")

# ---------------------------------------------------------------- números citados no texto
thf = np.radians(np.linspace(0.2, 180, 9000))
for abs_ in (False, True):
    ds = np.abs(nuclear_amplitude(res[abs_][0], 0.0, k, thf)) ** 2 * 10
    i = np.where((ds[1:-1] < ds[:-2]) & (ds[1:-1] < ds[2:]))[0] + 1
    print(("complexo" if abs_ else "real    "), "mínimos (graus, mb/sr):",
          [(round(float(np.degrees(thf[j])), 1), float(f"{ds[j]:.3g}")) for j in i],
          f"| dsig(0) = {ds[0]:.0f}, dsig(180) = {ds[-1]:.0f} mb/sr")
print("|S_l| (complexo), l = 0..12:", np.round(eta[:13], 3), "| mínimo em l =", lmin,
      "| b_l = (l+1/2)/k:", np.round((np.arange(13) + 0.5) / k, 2))
print(f"pi/(kR_V) = {np.degrees(np.pi/(k*RV)):.1f} graus; E_cm = {Ecm:.3f} MeV")
