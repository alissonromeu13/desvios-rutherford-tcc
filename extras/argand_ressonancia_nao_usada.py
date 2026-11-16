"""Figura NÃO usada no TCC (guardada como possível apêndice) – Diagrama de Argand e atraso de Wigner da ressonância de forma l = 3
(mesmo poço da Figura 12: K0R = 5,3; energias em hbar^2/(2 mu R^2), E = (kR)^2).

  T_l = e^{i delta_l} sin(delta_l) = (S_l - 1)/(2i)      -> círculo de unitariedade |T - i/2| = 1/2
  Breit-Wigner: delta = delta_bg + arctan[(Gamma/2)/(E_R - E)]
     T percorre o círculo no sentido anti-horário; em E_R +- Gamma/2 a fase ressonante vale pi/4 e 3pi/4
  Atraso de Wigner: Delta t = 2 hbar d delta/dE   (unidade: hbar/E0 = 2 mu R^2/hbar)
     BW: d delta/dE = (Gamma/2)/[(E - E_R)^2 + Gamma^2/4]  -> lorentziana de largura total Gamma
"""
import numpy as np
from scipy.special import spherical_jn as j, spherical_yn as y
from scipy.optimize import curve_fit
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from estilo import *
from matplotlib.ticker import FixedLocator, FixedFormatter

K0R = 5.3

def delta_mod(l, kR):
    KR = np.sqrt(kR**2 + K0R**2)
    num = kR * j(l, kR, derivative=True) * j(l, KR) - KR * j(l, kR) * j(l, KR, derivative=True)
    den = kR * y(l, kR, derivative=True) * j(l, KR) - KR * y(l, kR) * j(l, KR, derivative=True)
    return np.arctan2(num, den)

def delta_cont(l, kR, n_b):
    d = np.unwrap(2 * delta_mod(l, kR)) / 2
    return d - np.round(d[0] / np.pi) * np.pi + n_b * np.pi

# grade fina em energia
E = np.linspace(1e-4, 9.0, 90001); kR = np.sqrt(E)
d3 = delta_cont(3, kR, 0); d0 = delta_cont(0, kR, 2)
T3 = np.exp(1j * d3) * np.sin(d3)
S3 = np.exp(2j * d3)
assert np.allclose(np.abs(T3 - 0.5j), 0.5, atol=1e-12) and np.allclose(T3, (S3 - 1) / 2j)
print("unitariedade: |T_3 - i/2| = 1/2 em todos os pontos")

# ajuste de Breit-Wigner (mesmo procedimento da Figura 12)
bw = lambda E, ER, G, b0, b1: b0 + b1 * (E - 3.3) + np.arctan2(G / 2, ER - E)
m = (E > 2.9) & (E < 3.8)
(ER, G, b0, b1), _ = curve_fit(bw, E[m], d3[m], p0=(3.3, 0.4, 0, 0))
print(f"Breit-Wigner: E_R = {ER:.3f}, Gamma = {G:.3f}, fase de fundo em E_R = {b0:.3f} rad = {b0/np.pi:.3f} pi")

# atraso de Wigner
dt3 = 2 * np.gradient(d3, E); dt0 = 2 * np.gradient(d0, E)
ipk = np.argmax(dt3)
half = dt3[ipk] / 2
lo = E[:ipk][dt3[:ipk] < half][-1]; hi = E[ipk:][dt3[ipk:] < half][0]
print(f"atraso l=3: pico {dt3[ipk]:.2f} em E = {E[ipk]:.3f} (4/Gamma = {4/G:.2f}); largura a meia-altura = {hi-lo:.3f} (Gamma = {G:.3f})")
print(f"atraso l=0 no mesmo E: {dt0[ipk]:.3f} (negativo: avanço temporal, espalhamento não ressonante)")
assert abs((hi - lo) / G - 1) < 0.1

# pontos E_R +- Gamma/2 sobre o círculo
iR = np.argmin(np.abs(E - ER)); im_ = np.argmin(np.abs(E - (ER - G / 2))); ip_ = np.argmin(np.abs(E - (ER + G / 2)))
print(f"T_3 em E_R - Gamma/2: {T3[im_]:.3f}; em E_R: {T3[iR]:.3f}; em E_R + Gamma/2: {T3[ip_]:.3f}")
print(f"fração da volta percorrida em 0 < E < 9: {2*d3.max()/(2*np.pi):.2f}")

# ------------------------------------------------------------------ figura
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1, 1.25]})
th = np.linspace(0, 2 * np.pi, 400)
a1.plot(0.5 * np.cos(th), 0.5 + 0.5 * np.sin(th), color=CINZA, lw=0.8, ls="--")
a1.plot(T3.real, T3.imag, color=AMARELO, lw=1.2)
Ep = np.arange(0.1, 9.0, 0.1)                         # passos iguais de energia
ip = np.searchsorted(E, Ep)
a1.plot(T3.real[ip], T3.imag[ip], "o", ms=2.6, color=AMARELO, mec="none")
# atenção: com Re T = sin(2 delta)/2, o lado direito do círculo corresponde a delta = pi/4 (E < E_R)
marcas = [(im_, "$E_R - \\Gamma/2$", (7, 0), "left", "center"),
          (iR, "$E_R$", (0, 6), "center", "bottom"),
          (ip_, "$E_R + \\Gamma/2$", (-7, 0), "right", "center")]
for idx, rot, off, ha, va in marcas:
    a1.plot([T3.real[idx]], [T3.imag[idx]], "o", ms=5, color=TINTA, mec="white", mew=1, zorder=6)
    a1.annotate(rot, (T3.real[idx], T3.imag[idx]), xytext=off, textcoords="offset points",
                fontsize=7, ha=ha, va=va)
a1.plot([0], [0], "o", ms=4.5, mfc="white", mec=TINTA, mew=1, zorder=6)
a1.annotate("$E = 0$", (0, 0), xytext=(6, -2), textcoords="offset points", fontsize=7, va="top")
iE9 = len(E) - 1
a1.plot([T3.real[iE9]], [T3.imag[iE9]], "s", ms=4, mfc="white", mec=TINTA, mew=1, zorder=6)
a1.annotate("$kR = 3$", (T3.real[iE9], T3.imag[iE9]), xytext=(-6, -2), textcoords="offset points",
            fontsize=7, ha="right", va="top")
a1.annotate("", xy=(0.36, 0.86), xytext=(0.47, 0.70),
            arrowprops=dict(arrowstyle="->", lw=0.9, color=TINTA, connectionstyle="arc3,rad=0.25"))
a1.text(0.0, 0.5, "pontos a cada\n$\\Delta E = 0{,}1$", fontsize=6.8, ha="center", va="center", color=CINZA)
a1.set_aspect("equal")
a1.set(xlim=(-0.68, 0.68), ylim=(-0.12, 1.12), xlabel=r"$\mathrm{Re}\,T_3$", ylabel=r"$\mathrm{Im}\,T_3$")
a1.xaxis.set_major_locator(FixedLocator([-0.5, 0, 0.5])); a1.xaxis.set_major_formatter(FixedFormatter(["$-0{,}5$", "0", "0,5"]))
a1.yaxis.set_major_locator(FixedLocator([0, 0.5, 1])); a1.yaxis.set_major_formatter(FixedFormatter(["0", "0,5", "1"]))

# (b) atraso de Wigner
lor = 2 * (G / 2) / ((E - ER) ** 2 + G**2 / 4)
a2.axhline(0, color=CINZA, lw=0.6)
a2.plot(E, dt3, color=AMARELO, lw=1.5, ls=":", label="$l=3$ (exato)")
a2.plot(E, lor, color=TINTA, lw=0.9, ls="--", label="Breit–Wigner")
a2.plot(E, dt0, color=AZUL, lw=1.3, label="$l=0$ (exato)")
a2.annotate("", xy=(lo, half), xytext=(hi, half), arrowprops=dict(arrowstyle="<->", lw=0.8, color=TINTA, mutation_scale=7, shrinkA=0, shrinkB=0))
a2.text(hi + 0.05, half, "$\\Gamma \\approx " + f"{G:.2f}".replace(".", "{,}") + "$", fontsize=7.5, va="center")
a2.text(4.4, -1.3, "$l=0$: avanço ($\\Delta t < 0$)", fontsize=6.8, color=AZUL, ha="center", va="top")
a2.set(xlim=(1.5, 5.5), ylim=(-2.4, 11.5), xlabel=r"$E\ \,[\hbar^2/2\mu R^2]$",
       ylabel=r"$\Delta t = 2\hbar\,d\delta_l/dE\ \,[2\mu R^2/\hbar]$")
a2.legend(fontsize=7, loc="upper right")
virgula(a2)
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout()
salvar(fig, "extras/argand_ressonancia")
