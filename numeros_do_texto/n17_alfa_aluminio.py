"""Números da Seção 4.4.3 (alfa + 27Al, Figura 17) que não aparecem na própria figura.
Cada bloco imprime um número citado no texto, com a fonte dos parâmetros.
  1. Verificação da implementação: sigma_reação em 24,7 MeV para os três potenciais de McFadden e Satchler
     para o Al (1966, p. 187, Tab. 3), comparada aos valores publicados (1231, 1208 e 1159 mb).
  2. Robustez: chi2 dos seis pontos de Bieler (RaC e Po) com os potenciais do Al e com W reduzido.
  3. Funções de excitação em 90° e 135° (comparação com Rutherford e Chadwick, 1925, e Farwell e Wegner, 1954, p. 1212).
  4. Teste de absorção de superfície: a parte imaginária de volume é trocada por um termo de superfície
     (derivada de Woods-Saxon, forma de Koning e Delaroche 2003, Eq. 1, com pico W_D em R), mesma geometria.
     Não há parâmetros de superfície publicados para alfa + 27Al nessa energia; as profundidades são varridas.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from scipy.optimize import brentq
from scatter_core import *

A, Z = 27, 13
mT = A * AMU; mu = M_ALPHA * mT / (M_ALPHA + mT)
RC0 = 1.3                                  # raio coulombiano, r_C = 1,3 fm (McFadden e Satchler 1966, p. 179)
POT_TCC = (185.0, 25.0, 1.40, 0.52)        # V, W, r0, a (McFadden e Satchler 1966, p. 198)
POT_AL = [(49.1, 7.34, 1.701, 0.572, 1231), (127.8, 12.8, 1.534, 0.543, 1208), (197.1, 17.0, 1.349, 0.592, 1159)]  # Tab. 3, p. 187

def calc(Elab, V, W, r0, a, lmax=40, WD=0.0):
    Ecm = Elab * mT / (M_ALPHA + mT)
    k = np.sqrt(2 * mu * Ecm) / HBARC
    eta = 2 * Z * E2 * mu / (HBARC**2 * k)
    R = r0 * A ** (1 / 3); Rc = RC0 * A ** (1 / 3)
    sup = lambda r: 4 * np.exp((r - R) / a) / (1 + np.exp((r - R) / a)) ** 2      # -4a df/dr, vale 1 em r = R
    U = lambda r: coulomb_sphere(r, 2 * Z, Rc) - (V + 1j * W) * woods_saxon(r, R, a) - 1j * WD * sup(r)
    S, k = smatrix(U, mu, Ecm, lmax, max(R + 14 * a, 2.5 * eta / k + 2), eta=eta, h=0.01)
    return S, k, eta

def sigma_reacao(S, k):
    l = np.arange(len(S)); return np.pi / k**2 * np.sum((2 * l + 1) * (1 - np.abs(S) ** 2)) * 10   # mb

def razao(Elab, thl, pot=POT_TCC):
    S, k, eta = calc(Elab, *pot[:4], WD=pot[4] if len(pot) > 4 else 0.0)
    return ratio_to_rutherford(S, eta, k, lab_to_cm_angle(np.radians(np.atleast_1d(thl)), M_ALPHA, mT))

# diafragmas de Bieler (1924, p. 445, Tabs. I e II) e média com peso cossec^3(theta/2) (p. 438, Eq. 1; Apêndice F)
BINS = np.array([[17.6, 23.6], [23.7, 44.2], [43.9, 59.7], [60.3, 99.5]])
RAC = (np.array([1.00, 0.94, 0.78, 0.71]), np.array([0.03, 0.02, 0.02, 0.02]))
PO = (np.array([1.00, 0.93]), np.array([0.03, 0.03]))               # diafragmas III e IV
E_al = lambda R: 7.644 * (R / 7.0) ** (2 / 3)                       # Geiger, âncora do Apêndice F
E_RaC, E_Po = E_al(6.6), E_al(3.4)

def binned(Elab, pot):
    out = []
    for b in np.radians(BINS):
        tl = np.linspace(*b, 600); w = np.sin(tl / 2) ** -3
        r = razao(Elab, np.degrees(tl), pot)
        out.append(np.trapezoid(r * w, tl) / np.trapezoid(w, tl))
    return np.array(out)

def chi2_seis(pot):
    m1, m2 = binned(E_RaC, pot), binned(E_Po, pot)[2:]
    return np.sum(((m1 - RAC[0]) / RAC[1]) ** 2) + np.sum(((m2 - PO[0]) / PO[1]) ** 2), m1, m2

if __name__ == "__main__":
    print("1. sigma_reação em 24,7 MeV (potenciais do Al, Tab. 3):")
    for V, W, r0, a, pub in POT_AL:
        S, k, _ = calc(24.7, V, W, r0, a, lmax=60)
        s = sigma_reacao(S, k); print(f"   V = {V}: {s:.0f} mb (publicado {pub} mb; diferença {100*(s-pub)/pub:+.2f} %)")

    print(f"\n2. Bieler: E(RaC) = {E_RaC:.3f} MeV, E(Po) = {E_Po:.3f} MeV")
    c, m1, m2 = chi2_seis(POT_TCC)
    print(f"   potencial do TCC: RaC {np.round(m1, 3)}, Po III e IV {np.round(m2, 3)}, chi2 (6 pontos) = {c:.0f}")
    for E in (7.0, 7.68):
        print(f"   {E} MeV: diafragmas III e IV = {np.round(binned(E, POT_TCC)[2:], 2)}")
    for V, W, r0, a, _ in POT_AL:
        print(f"   potencial do Al V = {V}: chi2 = {chi2_seis((V, W, r0, a))[0]:.0f}")
    for W in (5.0, 10.0):
        c, m1, m2 = chi2_seis((185.0, W, 1.40, 0.52))
        print(f"   W = {W}: chi2 = {c:.0f}; RaC {np.round(m1, 2)}; Po {np.round(m2, 2)}")

    print("\n3. Funções de excitação:")
    for E in (5.0, 7.81):
        r = razao(E, [90.0, 135.0]); print(f"   {E} MeV: 90° -> {r[0]:.2f}, 135° -> {r[1]:.2f}")
    Es = np.arange(2.0, 6.0, 0.02); r135 = np.array([razao(E, 135.0)[0] for E in Es])
    print(f"   135°: sigma/sigma_R >= 0,98 até {Es[np.argmax(r135 < 0.98)]:.1f} MeV")
    # Rutherford e Chadwick (1925), alcances 4,0 / 5,0 / 5,3 / 6,7 / 6,8 cm
    for Rg in (4.0, 5.0, 5.3, 6.7, 6.8):
        E = E_al(Rg); r = razao(E, [90.0, 135.0])
        print(f"   alcance {Rg} cm (E = {E:.2f} MeV): 90° -> {r[0]:.2f}, 135° -> {r[1]:.2f}")
    print(f"   razão 135°: 5,3 cm / 6,8 cm = {razao(E_al(5.3), 135.0)[0] / razao(E_al(6.8), 135.0)[0]:.2f}")
    # Farwell e Wegner (1954, p. 1212): distâncias apsidais em 135° convertidas em energia e alcance
    thcm = lab_to_cm_angle(np.radians(135.0), M_ALPHA, mT); fac = 0.5 * (1 + 1 / np.sin(thcm / 2))
    for rmin in (13.0, 8.0, 7.0, 6.0):
        El = 2 * Z * E2 / (rmin / fac) * (M_ALPHA + mT) / mT
        print(f"   r_min = {rmin} fm: E_lab = {El:.2f} MeV, alcance = {7.0*(El/7.644)**1.5:.1f} cm, modelo 135° = {razao(El, 135.0)[0]:.2f}")

    print("\n4. Absorção puramente de superfície (V = 185 MeV, r0 = 1,40 fm, a = 0,52 fm, W de volume = 0):")
    for WD in (5.0, 10.0, 20.0, 40.0):
        pot = (185.0, 0.0, 1.40, 0.52, WD)
        c, m1, m2 = chi2_seis(pot)
        S, k, _ = calc(24.7, *pot[:4], lmax=60, WD=WD)
        print(f"   W_D = {WD:4.0f} MeV: chi2 = {c:.0f}; RaC {np.round(m1, 2)}; Po {np.round(m2, 2)};"
              f" sigma_reação(24,7 MeV) = {sigma_reacao(S, k):.0f} mb")
