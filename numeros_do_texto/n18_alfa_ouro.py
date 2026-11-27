"""Números da Seção 4.4.3 (alfa + 197Au, Figura 15) e do Apêndice C que não aparecem na própria figura.
Mesmo potencial da Figura 14 (McFadden e Satchler 1966, p. 198; r_C = 1,3 fm).
Dados de comparação: Farwell e Wegner (1954), Tab. I, p. 1216 (E_1/4 e K medidos no detector, ângulo +-1°).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from scipy.optimize import brentq
from scatter_core import *

A, Z = 197, 79
mT = A * AMU; mu = M_ALPHA * mT / (M_ALPHA + mT)
V, W, r0, a, rc = 185.0, 25.0, 1.40, 0.52, 1.3
R, Rc = r0 * A ** (1 / 3), rc * A ** (1 / 3)

def calc(Elab):
    Ecm = Elab * mT / (M_ALPHA + mT); k = np.sqrt(2 * mu * Ecm) / HBARC; eta = 2 * Z * E2 * mu / (HBARC**2 * k)
    U = lambda r: coulomb_sphere(r, 2 * Z, Rc) - (V + 1j * W) * woods_saxon(r, R, a)
    S, k = smatrix(U, mu, Ecm, int(k * (R + 6) + 15), max(R + 14 * a, 2.5 * eta / k + 2), eta=eta, h=0.01)
    return S, k, eta

def razao(E, thl):
    S, k, eta = calc(E)
    return ratio_to_rutherford(S, eta, k, lab_to_cm_angle(np.radians(np.atleast_1d(thl)), M_ALPHA, mT))

def fator_cin(thl, m=4.0026, M=196.97):
    """E_detector/E_incidente no espalhamento elástico (Apêndice C)."""
    t = np.radians(thl); return ((m * np.cos(t) + np.sqrt(M**2 - m**2 * np.sin(t) ** 2)) / (m + M)) ** 2

if __name__ == "__main__":
    th = np.linspace(20, 170, 151)
    print("Ausência de desvio (Seção 4.4.3):", {E: f"max|1 - razão| = {np.max(np.abs(1 - razao(E, th))):.1e}" for E in (4.7, 7.35, 7.68)})
    print("135°: 18 MeV ->", np.round(razao(18.0, 135.0), 2), "; 20 MeV ->", np.round(razao(20.0, 135.0), 2))
    print(f"r_min(60°)/r_min(180°) = {(1 + 1/np.sin(np.radians(30))) / 2:.2f}")
    Es = np.arange(12.0, 45.01, 0.25)
    r60 = np.array([razao(E, 60.0)[0] for E in Es])
    print(f"60°: primeira queda abaixo de 0,98 em {Es[np.argmax(r60 < 0.98)]:.1f} MeV; máximo {r60.max():.2f} em {Es[r60.argmax()]:.1f} MeV")
    print("Fatores cinemáticos (Apêndice C): 60° ->", f"{fator_cin(60):.3f}", "; 96° ->", f"{fator_cin(96):.3f}")
    # Farwell e Wegner (1954), Tab. I: E_1/4 (detector) e K (inclinação de ln sigma contra E)
    for thl, E14x, Kx in [(60, 33.6, 0.271), (96, 24.9, 0.394)]:
        f = fator_cin(thl)
        E14 = brentq(lambda E: razao(E, float(thl))[0] - 0.25, 20, 44)
        print(f"{thl}°: E_1/4 medido (incidente) = {E14x/f:.1f} +- {0.4/f:.1f} MeV; modelo = {E14:.1f} MeV")
    for thl in (59.0, 61.0):
        print(f"   ângulo efetivo {thl}°: modelo E_1/4 = {brentq(lambda E: razao(E, thl)[0] - 0.25, 30, 38):.1f} MeV")
    # inclinação K: sigma ∝ exp(-K E_det) acima de E_1/4; sigma = razão x sigma_R, sigma_R ∝ 1/E^2
    f = fator_cin(60); Es2 = np.arange(24.0, 45.01, 0.25)
    r = np.array([razao(E, 60.0)[0] for E in Es2]); Ec = f * Es2; lns = np.log(r) - 2 * np.log(Ec)
    for lo in (0.5, 0.25, 0.1):
        m = (r < lo) & (Ec <= 42.5)
        print(f"   60°, ajuste para razão < {lo} (E_det {Ec[m][0]:.1f}–{Ec[m][-1]:.1f} MeV): K = {-np.polyfit(Ec[m], lns[m], 1)[0]:.3f} /MeV (medido 0,271 +- 0,004)")
