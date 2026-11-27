"""Figura retirada da versão final do TCC (era a Figura 15, Seção 4.4.1); guardada em extras/. – Onda s no poço quadrado complexo U = -V0 - iW0, analítico.

 K = sqrt(2mu(E + V0 + iW0))/hbar,   S0 = e^{-2ikR} [K cot(KR) + ik] / [K cot(KR) - ik]
 sigma_el = (pi/k^2)|1 - S0|^2,  sigma_r = (pi/k^2)(1 - |S0|^2),  sigma_tot = (2pi/k^2)(1 - Re S0)  (teorema óptico)
 Limites: sigma_el, sigma_tot <= 4pi/k^2;  sigma_r <= pi/k^2.
 Probabilidade de absorção da onda s: 1 - eta0^2 = sigma_r k^2/pi  (máxima onde eta0 é mínimo).
Parâmetros ilustrativos tipo nêutron + 27Al: mu = 906 MeV, R = 3,6 fm, V0 = 40 MeV (E = E_cm).
 K0 R = 4,91 > 3pi/2: dois estados ligados s, o segundo raso -> a0 grande e positivo (como o dêuteron, Figura 10).
"""
import numpy as np
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from estilo import *
from scatter_core import HBARC
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

mu, R, V0 = 906.0, 3.6, 40.0
K0R = np.sqrt(2 * mu * V0) / HBARC * R
E = np.logspace(-2, np.log10(30), 3000)
k = np.sqrt(2 * mu * E) / HBARC

def calc(W0):
    K = np.sqrt(2 * mu * (E + V0 + 1j * W0) + 0j) / HBARC
    S0 = np.exp(-2j * k * R) * (K / np.tan(K * R) + 1j * k) / (K / np.tan(K * R) - 1j * k)
    sel = np.pi / k**2 * np.abs(1 - S0) ** 2 * 10          # mb
    sr = np.pi / k**2 * (1 - np.abs(S0) ** 2) * 10
    stot = 2 * np.pi / k**2 * (1 - S0.real) * 10
    return S0, sel, sr, stot

# ------------------------------------------------------------------ checagens
casos = [(0.0, AZUL, "-"), (1.0, LARANJA, "--"), (5.0, VERDE, "-."), (20.0, AMARELO, ":")]
res = {}
for W0, _, _ in casos:
    S0, sel, sr, stot = calc(W0)
    err = np.max(np.abs(stot - sel - sr) / stot); print(f"   W0 = {W0}: erro relativo máx. do teorema óptico = {err:.1e}"); assert err < 1e-8
    assert np.all(np.abs(S0) <= 1 + 1e-12)
    assert np.all(sel <= 4 * np.pi / k**2 * 10 + 1e-9) and np.all(sr <= np.pi / k**2 * 10 + 1e-9)
    res[W0] = (S0, sel, sr, stot)
print("teorema óptico sigma_tot = sigma_el + sigma_r (erro < 1e-8), |S0| <= 1 e limites unitários conferidos")
a0 = np.sqrt(res[0.0][1][0] / 10 / (4 * np.pi))
nb = int((K0R + np.pi / 2) // np.pi)
print(f"K0R = {K0R:.3f} (3pi/2 = {1.5*np.pi:.3f}): N_b = {nb} estados ligados s; a0 ~ +{a0:.2f} fm = {a0/R:.2f} R")
for W0 in (1.0, 5.0):
    sr = res[W0][2]; m = E < 0.05
    p = np.polyfit(np.log(E[m]), np.log(sr[m]), 1)[0]
    Eb = np.array([1e-6, 2e-6]); kb = np.sqrt(2 * mu * Eb) / HBARC
    Kb = np.sqrt(2 * mu * (Eb + V0 + 1j * W0) + 0j) / HBARC
    Sb = np.exp(-2j * kb * R) * (Kb / np.tan(Kb * R) + 1j * kb) / (Kb / np.tan(Kb * R) - 1j * kb)
    srb = (1 - np.abs(Sb) ** 2) / kb**2
    print(f"W0 = {W0}: inclinação de ln sigma_r x ln E: {np.log(srb[1]/srb[0])/np.log(2):.4f} em E ~ 1 eV, "
          f"{p:.3f} em E < 50 keV (lei 1/v assintótica: -0,5)")
minimos = {}
for W0 in (1.0, 5.0, 20.0):
    eta = np.abs(res[W0][0]); i = np.argmin(eta)
    minimos[W0] = (E[i], eta[i])
    print(f"W0 = {W0}: eta0 mínimo = {eta[i]:.3f} em E = {E[i]:.2f} MeV" + (" (borda do gráfico)" if i == len(E) - 1 else ""))
for W0 in (0.0, 1.0, 5.0, 20.0):
    sel = res[W0][1]; mm = (E > 5)
    print(f"W0 = {W0}: mínimo de sigma_el em E = {E[mm][np.argmin(sel[mm])]:.1f} MeV, valor {sel[mm].min():.1f} mb")

# identidade de fluxo (equação de continuidade, Seção 4.3.2): 1 - |S0|^2 = (8 mu / hbar^2 k) W0 int_0^R |u0|^2 dr
from scipy.integrate import quad
for W0 in (1.0, 5.0, 20.0):
    for Ei in (0.05, 2.1, 13.7):
        ki = np.sqrt(2 * mu * Ei) / HBARC; Ki = np.sqrt(2 * mu * (Ei + V0 + 1j * W0) + 0j) / HBARC
        Si = np.exp(-2j * ki * R) * (Ki / np.tan(Ki * R) + 1j * ki) / (Ki / np.tan(Ki * R) - 1j * ki)
        Ai = 0.5j * (np.exp(-1j * ki * R) - Si * np.exp(1j * ki * R)) / np.sin(Ki * R)
        Ii = quad(lambda r: abs(Ai * np.sin(Ki * r)) ** 2, 0, R, limit=200)[0]
        assert abs((1 - abs(Si) ** 2) - 8 * mu / (HBARC**2 * ki) * W0 * Ii) < 1e-9
print("identidade de fluxo 1 - |S0|^2 = (8 mu/hbar^2 k) W0 int |u0|^2 dr conferida")

# ------------------------------------------------------------------ figura
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.2, 3.15))
xt = FixedLocator([0.01, 0.1, 1, 10]); xf = FixedFormatter(["0,01", "0,1", "1", "10"])

# (a) teorema óptico, W0 = 5 MeV — legenda própria: aqui o estilo de linha distingue as grandezas
S0, sel, sr, stot = res[5.0]
h_lim, = a1.plot(E, 4 * np.pi / k**2 * 10, color=CINZA, lw=0.8, label=r"$4\pi/k^2$ (limite unitário)")
h_tot, = a1.plot(E, stot, color=TINTA, lw=1.8, label=r"$\sigma_0^{\mathrm{tot}} = \sigma_0^{\mathrm{elás}} + \sigma_0^{\mathrm{reação}}$")
h_el, = a1.plot(E, sel, color=VERDE, lw=1.4, ls="-", label=r"$\sigma_0^{\mathrm{elás}}$")
h_r, = a1.plot(E, sr, color=VERDE, lw=1.4, ls=(0, (2.2, 1.6)), label=r"$\sigma_0^{\mathrm{reação}}$")
idx = np.searchsorted(E, np.logspace(-1.7, np.log10(25), 9))
h_ot, = a1.plot(E[idx], (2 * np.pi / k[idx]**2 * (1 - S0[idx].real)) * 10, "o", ms=3.8, mfc="white", mec=TINTA,
                mew=0.9, zorder=5, ls="none", label=r"$(2\pi/k^2)(1-\mathrm{Re}\,S_0)$")
a1.legend(handles=[h_tot, h_ot, h_el, h_r, h_lim], loc="upper right", fontsize=6.2, title="$W_0 = 5$ MeV",
          title_fontsize=6.8, handlelength=2.2, borderaxespad=0.3, labelspacing=0.3, frameon=True, framealpha=0.92,
          edgecolor="none")
a1.text(0.04, 0.04, r"lei $1/v$ ($E \to 0$):" + "\n" + r"$\sigma_0^{\mathrm{reação}} \propto E^{-1/2}$", transform=a1.transAxes,
        fontsize=6.8, color=CINZA)
a1.set(xscale="log", yscale="log", ylim=(10, 3e7), xlim=(0.01, 30), xlabel=r"$E_\mathrm{cm}$ [MeV]", ylabel=r"$\sigma_0$ [mb]")

# (b) sigma_el para vários W0
handles = []
for W0, c, ls in casos:
    h, = a2.plot(E, res[W0][1], color=c, ls=ls, lw=1.4, label=f"$W_0 = {W0:g}$ MeV")
    handles.append(h)
a2.annotate("estado ligado raso\n" + (f"$a_0 \\approx +{a0:.1f}$ fm").replace(".", "{,}"), xy=(0.012, res[0.0][1][1]),
            xytext=(0.014, 90), fontsize=6.8, arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA))
a2.text(13, 2.5e3, "mínimo\nelástico", fontsize=6.8, ha="center", va="bottom")
a2.annotate("", xy=(13.0, 60), xytext=(13.0, 2.4e3), arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA))
a2.set(xscale="log", yscale="log", ylim=(10, 2e4), xlim=(0.01, 30), xlabel=r"$E_\mathrm{cm}$ [MeV]",
       ylabel=r"$\sigma_0^{\mathrm{elás}}$ [mb]")

# (c) eta0
for W0, c, ls in casos:
    a3.plot(E, np.abs(res[W0][0]), color=c, ls=ls, lw=1.4)
for W0, c in [(1.0, LARANJA), (5.0, VERDE)]:
    Em, em = minimos[W0]
    a3.plot([Em], [em], "o", ms=5, color=c, mec="white", mew=1, zorder=5)
a3.annotate("absorção\nmáxima", xy=minimos[5.0], xytext=(0.03, 0.42), fontsize=6.8,
            arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA, connectionstyle="arc3,rad=0.2"))
a3.annotate("", xy=minimos[1.0], xytext=(0.12, 0.55), arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA,
            connectionstyle="arc3,rad=-0.25"))
a3.text(29, 0.33, r"$W_0 = 20$: mínimo" + "\nacima de 30 MeV", fontsize=6.5, ha="right", va="top", color=AMARELO)
a3.axhline(1, color=CINZA, lw=0.7, ls="--")
a3.set(xscale="log", ylim=(0, 1.05), xlim=(0.01, 30), xlabel=r"$E_\mathrm{cm}$ [MeV]", ylabel=r"$\eta_0 = |S_0|$")
virgula(a3, "y")
for ax in (a1, a2, a3):
    ax.xaxis.set_major_locator(xt); ax.xaxis.set_major_formatter(xf); ax.xaxis.set_minor_formatter(NullFormatter())
rotulo(a1, "a"); rotulo(a2, "b"); rotulo(a3, "c")
fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=0.8)
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
salvar(fig, "extras/fig06_poco_complexo_absorcao_v4")
