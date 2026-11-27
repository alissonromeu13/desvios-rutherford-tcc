"""Figura 10 do TCC (Seção 4.2.2) – Comprimento de espalhamento da onda s no poço quadrado esférico.
 Painel único (v1.2): a seção de choque em função de kR foi retirada por exceder o regime kR << 1 da seção.

 delta_0 = -kR + arctan[(k/K) tan(KR)],  K = sqrt(k^2 + K0^2),  K0^2 = 2 mu V0/hbar^2
 a_0 = R [1 - tan(K0 R)/(K0 R)]                     (comprimento de espalhamento)
 r_e = R - R^3/(3 a_0^2) - 1/(K0^2 a_0)              (alcance efetivo do poço quadrado)
 k cot(delta_0) ~ -1/a_0 + r_e k^2/2                 (expansão de alcance efetivo, ERE)
 sigma_0 = (4 pi / k^2) sin^2 delta_0 <= 4 pi/k^2;   k->0: 4 pi a_0^2;  esfera rígida: a_0 = R
 Estados ligados de onda s: N_b = número de valores (n+1/2) pi < K0 R;  Levinson: delta_0(0) = N_b pi.
 Sistema n-p: poço quadrado equivalente a cada canal, ajustado a (a, r_e) experimentais
 (singleto: a = -23,74 fm, r = 2,77 fm; tripleto: a = +5,419 fm, r = 1,753 fm; Babenko e Petrov 2007, p. 672, Eq. 17,
 que reproduzem a compilação de Machleidt 2001; mu = m_N/2).
"""
import numpy as np
from scipy.optimize import brentq
from estilo import *
from scatter_core import HBARC

a0_de = lambda x: 1 - np.tan(x) / x                       # a0/R em função de x = K0 R
re_de = lambda x: 1 - 1 / (3 * a0_de(x) ** 2) - 1 / (x**2 * a0_de(x))   # r_e/R

def delta0(kR, K0R):
    KR = np.sqrt(kR**2 + K0R**2)
    return -kR + np.arctan(kR / KR * np.tan(KR))

# ------------------------------------------------------------------ checagens
for K0R in (1.0, 1.45, 2.0, 4.0):
    k = np.array([1e-3, 2e-3]); y = k / np.tan(delta0(k, K0R))
    r_num = 2 * (y[1] - y[0]) / (k[1] ** 2 - k[0] ** 2)
    assert abs(r_num - re_de(K0R)) < 1e-4, K0R
    kk = np.linspace(1e-3, 3e-3, 2); assert np.allclose(kk / np.tan(delta0(kk, K0R)), -1 / a0_de(K0R), atol=1e-4)
print("alcance efetivo: fórmula fechada confere com o ajuste numérico de k cot(delta0)")
for K0R in (1.0, 1.45, 2.0, 4.0, 5.3, 8.0):
    kR = np.linspace(1e-4, 80, 400000)
    d = np.unwrap(2 * delta0(kR, K0R)) / 2
    d0 = d[0] - np.round(d[-1] / np.pi) * np.pi
    Nb = int((K0R + np.pi / 2) // np.pi)
    assert abs(d0 / np.pi - Nb) < 1e-3
print("Levinson: delta0(0) = N_b pi conferido para K0R = 1; 1,45; 2; 4; 5,3; 8")

MU_NP = 938.919 / 2
def poco_equivalente(a, r, lo, hi):
    """(a, r_e) experimentais -> (K0 R, R) do poço quadrado no ramo [lo, hi]."""
    R_de = lambda x: a / a0_de(x)
    g = lambda x: R_de(x) * re_de(x) - r
    xs = np.linspace(lo, hi, 4000); gs = g(xs)
    i = np.where(np.sign(gs[:-1]) != np.sign(gs[1:]))[0][0]
    x = brentq(g, xs[i], xs[i + 1]); return x, R_de(x)
np_canais = {}
for nome, a, r, lo, hi in [("singleto", -23.74, 2.77, 1.0, np.pi / 2 - 1e-6),
                           ("tripleto", 5.419, 1.753, np.pi / 2 + 1e-6, 2.5)]:
    x, R = poco_equivalente(a, r, lo, hi)
    K0 = x / R; V0 = (K0 * HBARC) ** 2 / (2 * MU_NP)
    s = f"n-p {nome}: K0R = {x:.4f}, R = {R:.3f} fm, V0 = {V0:.1f} MeV, a/R = {a/R:+.3f}"
    if a > 0:
        f = lambda kap: np.sqrt(K0**2 - kap**2) / np.tan(np.sqrt(K0**2 - kap**2) * R) + kap
        kap = brentq(f, 1e-4, 0.999 * K0); EB = (kap * HBARC) ** 2 / (2 * MU_NP)
        s += f", E_B previsto = {EB:.3f} MeV (experimental 2,2246 MeV)"
    print(s)
    np_canais[nome] = (x, a / R)

# ------------------------------------------------------------------ figura
fig, a1 = plt.subplots(figsize=(4.6, 3.4))
ymin, ymax = -11, 8
polos = [(n + 0.5) * np.pi for n in range(3)]
limites = [0] + polos + [3 * np.pi]
for n in range(4):                                         # faixas de N_b
    if n % 2 == 1:
        a1.axvspan(limites[n], limites[n + 1], color=CINZA, alpha=0.07, lw=0)
    a1.text(0.5 * (limites[n] + limites[n + 1]), ymax - 0.3, f"$N_b = {n}$", ha="center", va="top", fontsize=7.5, color=TINTA)
for xp in polos:
    a1.axvline(xp, color=LARANJA, lw=0.8, ls=":")
for n in range(4):                                         # ramos separados: sem traço vertical no polo
    x = np.linspace(limites[n] + 1e-4, limites[n + 1] - 1e-4, 3000)
    a1.plot(x, a0_de(x), color=AZUL, lw=1.5)
a1.axhline(0, color=TINTA, lw=0.6)
a1.axhline(1, color=CINZA, ls="--", lw=1)
a1.text(5.2, 2.0, "$a_0 = R$", fontsize=7, color=CINZA, va="bottom")
for xc in (np.pi, 2 * np.pi, 3 * np.pi):                  # cruzamentos a0 = R
    a1.plot([xc], [1], "o", ms=4.5, mfc="white", mec=CINZA, mew=1, zorder=4, clip_on=False)
for K0R, c in [(1.0, AZUL), (1.45, LARANJA), (np.pi, VERDE)]:
    a1.plot([K0R], [a0_de(K0R)], "o", color=c, ms=5, mec="white", zorder=5)
# estados no limiar
a1.annotate("estado virtual", xy=(np.pi / 2 - 0.04, -7.2), xytext=(2.0, -6.2),
            fontsize=7, color=TINTA, va="center", arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA))
a1.annotate("estado ligado raso", xy=(np.pi / 2 + 0.06, 6.0), xytext=(2.1, 4.6),
            fontsize=7, color=TINTA, va="center", arrowprops=dict(arrowstyle="->", lw=0.6, color=TINTA))
# sistema nêutron-próton
xs_, ys_ = np_canais["singleto"]; xt_, yt_ = np_canais["tripleto"]
a1.plot([xs_], [ys_], "D", ms=5, color=TINTA, mec="white", zorder=6)
a1.plot([xt_], [yt_], "s", ms=5, color=TINTA, mec="white", zorder=6)
a1.text(xs_ + 0.28, ys_, r"$^1S_0$", fontsize=7, va="center",
        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
a1.text(xt_ + 0.15, yt_ + 0.25, r"$^3S_1$ (dêuteron)", fontsize=7, va="bottom",
        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
a1.set(xlabel=r"$K_0R=\sqrt{2\mu V_0}\,R/\hbar$", ylabel=r"$a_0/R$", ylim=(ymin, ymax), xlim=(0, 3 * np.pi))
a1.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi, 5 * np.pi / 2, 3 * np.pi])
a1.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$", r"$5\pi/2$", r"$3\pi$"])
from matplotlib.ticker import FixedLocator
a1.yaxis.set_major_locator(FixedLocator([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8]))
virgula(a1, "y")

# checagem numérica da expansão de alcance efetivo (não desenhada)
kR = np.linspace(1e-3, 4, 1600)
K0R = 1.45; a, r = a0_de(K0R), re_de(K0R)
kcot = -1 / a + 0.5 * r * kR**2
s_ere = 4 / (kR**2 + kcot**2)
s_ex = 4 * np.sin(delta0(kR, K0R)) ** 2 / kR**2
dev = np.abs(s_ere / s_ex - 1)
k5 = kR[np.argmax(dev > 0.05)]
print(f"K0R=1,45: a0/R = {a:.3f}, r_e/R = {r:.4f}; ERE desvia mais de 5% do exato a partir de kR = {k5:.2f}")
fig.tight_layout()
salvar(fig, "6.3/fig03_onda_s_comprimento_espalhamento_v2")
for K0R in [1.0, 1.45]:
    k = 1e-4
    print(f"K0R={K0R}: sigma0(k->0)/piR^2 = {4*np.sin(delta0(k, K0R))**2/k**2:.4f}  vs 4 a0^2 = {4*a0_de(K0R)**2:.4f}")
