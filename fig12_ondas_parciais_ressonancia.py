"""Figura 11 do TCC (Seção 4.2.3) – Seções de choque parciais e ressonância de forma. Poço quadrado atrativo,
K0R = 5,3 (analítico). Unidades: kR adimensional; energias em hbar^2/(2 mu R^2), isto é, E -> (kR)^2.

  tan delta_l = [k j_l'(kR) j_l(KR) - K j_l(kR) j_l'(KR)] / [k y_l'(kR) j_l(KR) - K y_l(kR) j_l'(KR)]
  sigma_l = (4 pi/k^2)(2l+1) sin^2 delta_l  <=  (4 pi/k^2)(2l+1)   (limite de unitariedade)
  Levinson: delta_l(0) = N_b(l) pi;  N_b(0) = # {(n+1/2)pi < K0R},  N_b(l>=1) = # zeros de j_{l-1} abaixo de K0R
  Breit-Wigner: delta_3 = delta_bg + arctan[(Gamma/2)/(E_R - E)];  atraso de Wigner Delta t = 2 hbar d delta/dE
  Potencial efetivo (inset): (2 mu R^2/hbar^2) V_ef = -(K0R)^2 Theta(1 - r/R) + l(l+1)/(r/R)^2
"""
import numpy as np
from scipy.special import spherical_jn as j, spherical_yn as y
from scipy.optimize import curve_fit
from estilo import *
from matplotlib.ticker import FixedLocator, FixedFormatter

K0R = 5.3
LMAX = 5

def delta_mod(l, kR):
    """delta_l definido a menos de múltiplos de pi."""
    KR = np.sqrt(kR**2 + K0R**2)
    num = kR * j(l, kR, derivative=True) * j(l, KR) - KR * j(l, kR) * j(l, KR, derivative=True)
    den = kR * y(l, kR, derivative=True) * j(l, KR) - KR * y(l, kR) * j(l, KR, derivative=True)
    return np.arctan2(num, den)

def n_ligados(l):
    if l == 0:
        return int((K0R + np.pi / 2) // np.pi)
    x = np.linspace(1e-3, K0R, 400000); f = j(l - 1, x)
    return int(np.sum(np.sign(f[:-1]) != np.sign(f[1:])))

def delta_abs(l, kR):
    """Fase contínua com a normalização de Levinson, delta_l(0) = N_b(l) pi."""
    d = np.unwrap(2 * delta_mod(l, kR)) / 2
    return d - np.round(d[0] / np.pi) * np.pi + n_ligados(l) * np.pi

# ------------------------------------------------------------------ checagens
for l in range(1, 5):                                     # lei de limiar de Wigner
    k = np.array([1e-2, 2e-2]); t = np.abs(np.tan(delta_mod(l, k)))
    assert abs(np.log(t[1] / t[0]) / np.log(2) - (2 * l + 1)) < 0.01
print("lei de Wigner tan(delta_l) ~ k^(2l+1) conferida para l = 1..4")
kk = np.linspace(1e-3, 120, 600000)
for l in range(LMAX + 1):                                  # Levinson: delta_l(infinito) -> 0
    d = delta_abs(l, kk)
    assert abs(d[-1]) < 0.15, (l, d[-1])
print("Levinson conferido (delta_l(0) = N_b pi e delta_l -> 0 em k grande):",
      {l: n_ligados(l) for l in range(LMAX + 1)})
x = np.linspace(0.1, 10, 200000); f = j(2, x)
K0R_lim = x[np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]]
print(f"o nível l=3 se tornaria ligado em K0R = {K0R_lim:.3f} (primeiro zero de j_2); aqui K0R = {K0R}")

kR = np.linspace(0.02, 3.0, 6000)
deltas = {l: delta_abs(l, kR) for l in range(LMAX + 1)}
sig = {l: 4 * (2 * l + 1) * np.sin(deltas[l]) ** 2 / kR**2 for l in deltas}   # em unidades de pi R^2
for l in sig:
    assert np.all(sig[l] <= 4 * (2 * l + 1) / kR**2 + 1e-12)
print("limite de unitariedade respeitado para l = 0..5")
tot = sum(sig.values())

# ressonância l = 3: pico, cruzamento de pi/2, ajuste Breit-Wigner, atraso de Wigner
kf = np.linspace(1.5, 2.2, 40001); d3 = delta_abs(3, np.concatenate([[0.02], kf]))[1:]
s3 = 28 * np.sin(d3) ** 2 / kf**2
ip = np.argmax(s3); k_pico = kf[ip]
k_cruz = kf[np.where((d3[:-1] < np.pi / 2) & (d3[1:] >= np.pi / 2))[0][0]]
E = kf**2
bw = lambda E, ER, G, b0, b1: b0 + b1 * (E - 3.3) + np.arctan2(G / 2, ER - E)
m = (kf > 1.70) & (kf < 1.95)
(ER, G, b0, b1), _ = curve_fit(bw, E[m], d3[m], p0=(3.33, 0.4, 0.0, 0.0))
dd = np.gradient(d3, E); dt_max = 2 * dd.max()                 # Delta t em unidades de 2 mu R^2/hbar... ver abaixo
t_livre = 1 / (2 * k_cruz)                                      # R/v = mu R/(hbar k) na mesma unidade
print(f"ressonância l=3: pico de sigma_3 em kR = {k_pico:.4f} (sigma_3 = {s3[ip]:.3f}, limite {28/k_pico**2:.3f}); "
      f"delta_3 = pi/2 em kR = {k_cruz:.4f}")
print(f"Breit-Wigner: E_R = {ER:.3f}, Gamma = {G:.3f} (unid. hbar^2/2muR^2), Gamma/E_R = {G/ER:.3f}; "
      f"barreira em r=R: l(l+1) = 12 (kR = {np.sqrt(12):.2f})")
print(f"fase de fundo no ajuste de Breit-Wigner à fase: b0 = {b0:.2f} rad (Seção 4.2.3.1)")
print(f"atraso de Wigner no pico: Delta t = {dt_max:.2f} vs passagem livre R/v = {t_livre:.3f}  ->  razão {dt_max/t_livre:.0f}")
assert abs(4 / G - dt_max) / dt_max < 0.05                      # BW: d delta/dE máx = 2/Gamma
meia = E[s3 >= s3[ip] / 2]                                      # largura a meia-altura de sigma_3(E), citada na Seção 4.2.3.1
print(f"largura a meia-altura de sigma_3: Gamma = {meia.max() - meia.min():.3f}; E no pico = {k_pico**2:.3f}; "
      f"Gamma/E_R = {(meia.max() - meia.min()) / k_pico**2:.3f}")

# ------------------------------------------------------------------ figura
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.1, 1.1, 0.85]})
cores = [AZUL, LARANJA, VERDE, AMARELO, CINZA]
handles = []
for l in range(LMAX):
    h, = a1.plot(kR, sig[l], color=cores[l], ls=LS[l % 4], lw=1.3, label=f"$l={l}$")
    a2.plot(kR, deltas[l] / np.pi, color=cores[l], ls=LS[l % 4], lw=1.3)
    handles.append(h)
h, = a1.plot(kR, tot, color=TINTA, lw=1.0, label=r"total ($l\leq5$)"); handles.append(h)

# (a) seções de choque
a1.plot(kR, 28 / kR**2, color=AMARELO, lw=0.8, ls=(0, (1, 2)))
a1.text(2.25, 6.4, "limite de\nunitariedade ($l=3$)", fontsize=6.8, color=CINZA, va="bottom")
a1.plot([k_pico], [s3[ip]], "o", ms=5, color=TINTA, mec="white", mew=1, zorder=6)
txt = ("ressonância $l=3$\n"
       + f"$kR_R \\approx {k_pico:.2f}$".replace(".", "{,}"))
a1.annotate(txt, xy=(k_pico, s3[ip]), xytext=(0.15, 12.8), fontsize=7, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=CINZA, lw=0.5),
            arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA, connectionstyle="arc3,rad=-0.2"))
a1.set(xlabel=r"$kR$", ylabel=r"$\sigma_l/\pi R^2$", ylim=(0, 14), xlim=(0, 3))

# (b) fases absolutas (Levinson)
for yv in (0.5, 1.5):
    a2.axhline(yv, color=CINZA, lw=0.7, ls="--")
a2.text(2.97, 0.47, r"$\mathrm{sen}^2\delta_l = 1$", fontsize=6.8, color=CINZA, ha="right", va="top")
for yv, t in [(2, r"$N_b = 2$ ($l=0$)"), (1, r"$N_b = 1$ ($l=1, 2$)"), (0, r"$N_b = 0$ ($l\geq3$)")]:
    a2.text(0.06, yv + 0.05, t, fontsize=6.8, va="bottom")
a2.plot([k_cruz], [0.5], "o", ms=5, color=TINTA, mec="white", mew=1, zorder=6)
a2.text(k_cruz - 0.1, 0.56, "$\\delta_3$ sobe:\nressonância", fontsize=6.8, ha="right", va="bottom")
i0 = np.argmin(np.abs(deltas[0] - 1.5 * np.pi))
a2.plot([kR[i0]], [1.5], "o", ms=4.5, mfc="white", mec=AZUL, mew=1, zorder=6)
a2.text(kR[i0] + 0.1, 1.6, "$\\delta_0$ desce:\nsem ressonância", fontsize=6.8, va="bottom")
a2.set(xlabel=r"$kR$", ylabel=r"$\delta_l/\pi$", xlim=(0, 3), ylim=(-0.08, 2.3))
a2.yaxis.set_major_locator(FixedLocator([0, 0.5, 1, 1.5, 2]))
a2.yaxis.set_major_formatter(FixedFormatter(["0", "0,5", "1", "1,5", "2"]))

# (c) potencial efetivo l = 3
xr = np.linspace(0.4, 3.5, 1500)
U = np.where(xr < 1, -K0R**2, 0.0) + 12 / xr**2
x_in, x_out = np.sqrt(12 / (ER + K0R**2)), np.sqrt(12 / ER)
a3.axvspan(1, x_out, color=CINZA, alpha=0.13, lw=0)
a3.axhline(0, color=CINZA, lw=0.6)
a3.plot(xr, U, color=AMARELO, lw=1.5)
a3.plot([x_in, x_out], [ER, ER], color=TINTA, lw=1.0, ls="--")
a3.plot([1], [12], "o", ms=4, color=AMARELO, mec="white", mew=0.8, zorder=5)
a3.text(1.08, 12.3, "barreira: $l(l+1) = 12$", fontsize=6.8, va="bottom")
a3.text(x_out + 0.08, ER + 0.4, "$E_R$", fontsize=7.5, va="bottom")
a3.text(0.5 * (1 + x_out), -5.5, "túnel", fontsize=6.8, color=CINZA, ha="center")
a3.text(0.46, -27.5, "poço", fontsize=6.8, color=CINZA)
a3.text(3.4, -27.5, "$l=3$", fontsize=7.5, ha="right")
a3.set(xlim=(0.4, 3.5), ylim=(-31, 16), xlabel=r"$r/R$", ylabel=r"$V_\mathrm{efetivo}\ \,[\hbar^2/2\mu R^2]$")
a3.set_xticks([1, 2, 3])

virgula(a1); virgula(a2, "x"); virgula(a3, "x")
rotulo(a1, "a"); rotulo(a2, "b"); rotulo(a3, "c")
fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=1.0)
fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
salvar(fig, "6.3/fig04_ondas_parciais_ressonancia_v3")
print(f"pontos de retorno em E_R: r/R = {x_in:.3f} (interno) e {x_out:.3f} (externo); vida média hbar/Gamma = {1/G:.2f} "
      f"= {1/G/t_livre:.1f} R/v")
print(f"delta_0 passa por 3pi/2 (descendo) em kR = {kR[i0]:.3f}")
