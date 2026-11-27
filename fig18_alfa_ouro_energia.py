"""Figura 15 do TCC (Seção 4.4.3) – Por que Rutherford funciona no ouro e quando deixa de funcionar.
alfa + 197Au com o mesmo potencial óptico (McFadden & Satchler, 1966; R_C = 1,3 A^{1/3} fm).
(a) sigma/sigma_R x theta para várias energias; (b) sigma/sigma_R em theta_lab = 60°, 96° e 135° em função da energia,
análogo às 'funções de excitação' de Farwell & Wegner (1954), com os E_1/4 medidos por eles (Tab. I, p. 1216).
"""
import numpy as np
from estilo import *
from scatter_core import *

A, Z = 197, 79
mT = A * AMU; mu = M_ALPHA * mT / (M_ALPHA + mT)
V, W, r0, a, rc = 185.0, 25.0, 1.40, 0.52, 1.3

def calc(Elab, lmax=None):
    Ecm = Elab * mT / (M_ALPHA + mT)
    k = np.sqrt(2 * mu * Ecm) / HBARC
    eta = 2 * Z * E2 * mu / (HBARC**2 * k)
    R = r0 * A ** (1 / 3); Rc = rc * A ** (1 / 3)
    lmax = lmax or int(k * (R + 6) + 15)
    U = lambda r: coulomb_sphere(r, 2 * Z, Rc) - (V + 1j * W) * woods_saxon(r, R, a)
    rm = max(R + 14 * a, 2.5 * eta / k + 2)
    S, k = smatrix(U, mu, Ecm, lmax, rm, eta=eta, h=0.01)
    return S, k, eta

def ratio_lab(S, k, eta, thl):
    return ratio_to_rutherford(S, eta, k, lab_to_cm_angle(np.radians(thl), M_ALPHA, mT))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
thl = np.linspace(5, 175, 400)
for E, c, ls in [(5.0, CINZA, ":"), (7.68, AZUL, "-"), (22.0, LARANJA, "--"), (30.0, VERDE, "-."), (40.0, AMARELO, (0, (5, 1.5, 1, 1.5)))]:
    S, k, eta = calc(E)
    S2, k2, e2 = calc(E, lmax=len(S) + 20)
    dif = np.max(np.abs(ratio_lab(S, k, eta, thl) - ratio_lab(S2, k2, e2, thl)))
    r = ratio_lab(S, k, eta, thl)
    print(f"E={E}: eta={eta:.1f} lmax={len(S)-1} conv={dif:.1e}  sigma/sigmaR(175°)={r[-1]:.3f}")
    a1.plot(thl, r, color=c, ls=ls, zorder=5 if E == 5.0 else 2, lw=2.0 if E == 5.0 else 1.6, label=(f"{E:g} MeV").replace(".", ","))
a1.set(yscale="log", ylim=(1e-3, 3), xlim=(0, 180), xlabel=r"$\theta_\mathrm{lab}$ [graus]", ylabel=r"$\sigma/\sigma_\mathrm{R}$")
a1.set_xticks([0, 30, 60, 90, 120, 150, 180])
a1.legend(fontsize=7.5, loc="lower left", title=r"$\alpha$ + $^{197}$Au", title_fontsize=8)

Es = np.linspace(10, 45, 71)
for thf, ls in [(135, "-"), (96, (0, (5, 1.5))), (60, "-.")]:
    r = [ratio_lab(*calc(E), np.array([float(thf)]))[0] for E in Es]
    a2.plot(Es, r, color=TINTA, ls=ls)
# rótulos diretos nas curvas (a legenda ficava sobre a curva de 60°)
for txt, x, y in [(r"$135^\circ$", 29.6, 0.019), (r"$96^\circ$", 36.4, 0.019), (r"$60^\circ$", 38.6, 0.16)]:
    a2.text(x, y, txt, fontsize=8, color=TINTA, ha="left", va="center")
# Farwell e Wegner (1954), Tab. I (p. 1216): E_1/4 para Au em 60° e 96°, energias medidas no detector,
# ~2 % abaixo da incidente por causa do recuo (p. 1215). Convertidas para energia incidente pelo fator cinemático.
def fator_cin(thl, m=4.0026, M=196.97):
    t = np.radians(thl)
    return ((m * np.cos(t) + np.sqrt(M**2 - m**2 * np.sin(t)**2)) / (m + M))**2
fw = [(60, 33.6, 0.4), (96, 24.9, 0.4)]
Efw = np.array([E / fator_cin(t) for t, E, _ in fw]); dEfw = np.array([d / fator_cin(t) for t, _, d in fw])
print("Farwell-Wegner E_1/4 incidente:", np.round(Efw, 2), "+-", np.round(dEfw, 2))
a2.errorbar(Efw, [0.25, 0.25], xerr=dEfw, fmt="o", color=TINTA, mfc="white", ms=5, capsize=2, lw=0.9, zorder=6,
            label="Farwell e Wegner (1954)")
# barreira no raio de toque (estimativa) e barreira l = 0 do próprio potencial óptico, como na Figura 14
Eb = 2 * Z * E2 / (1.2 * (4 ** (1 / 3) + A ** (1 / 3))) * (M_ALPHA + mT) / mT
from scipy.optimize import minimize_scalar
R0, Rc0 = r0 * A ** (1 / 3), rc * A ** (1 / 3)
Vef = lambda r: coulomb_sphere(r, 2 * Z, Rc0) - V * woods_saxon(r, R0, a)
opt = minimize_scalar(lambda r: -Vef(r), bounds=(7, 16), method="bounded")
EbO = Vef(opt.x) * (M_ALPHA + mT) / mT
print(f"barreira de toque: E_lab = {Eb:.2f} MeV")
print(f"barreira l=0 do potencial óptico: r_B = {opt.x:.2f} fm, V_B(cm) = {Vef(opt.x):.2f} MeV, E_lab = {EbO:.2f} MeV")
for Ex, ls, yt, txt in [(EbO, "--", 0.06, "barreira $l=0$\ndo potencial\nóptico:\n"), (Eb, ":", 0.0115, "barreira\n(toque):\n")]:
    a2.axvline(Ex, color=CINZA, lw=0.8, ls=ls)
    a2.text(Ex - 0.4, yt, txt + (f"${Ex:.1f}$ MeV").replace(".", "{,}"), fontsize=7, color=CINZA, ha="right", va="bottom")
a2.axhline(1, color=CINZA, lw=0.8, ls=":")
a2.set(yscale="log", ylim=(1e-2, 2), xlim=(10, 45), xlabel=r"$E_\alpha$ (lab) [MeV]", ylabel=r"$\sigma/\sigma_\mathrm{R}$")
a2.legend(fontsize=7, loc="upper right", handlelength=1.2, borderaxespad=0.3)
virgula(a1, "x"); virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout()
salvar(fig, "6.4/fig09_alfa_au_dependencia_energia")
