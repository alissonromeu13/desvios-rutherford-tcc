"""Fator de forma de carga confrontado com Bieler (Seção 4.1.2).
Na primeira aproximação de Born, sigma/sigma_R = |F(q)|^2. Esfera homogênea, R = 1,2 A^(1/3) fm,
partículas do RaC no alumínio com E_lab = 7,35 MeV (Apêndice F), k = 1,033 fm^-1 no centro de massa,
ângulos efetivos de laboratório dos quatro diafragmas (Tabela do Apêndice F) convertidos ao centro de massa (Apêndice E)."""
import numpy as np

tau = 4.0026 / 26.9815          # m_alfa / M_Al
k = 1.033                        # fm^-1 (E_cm = 6,40 MeV)
R = 1.2 * 27 ** (1 / 3)          # fm
medido = [1.00, 0.94, 0.78, 0.71]
for th_lab, r in zip((20.2, 30.9, 50.7, 75.8), medido):
    t = np.radians(th_lab)
    th_cm = t + np.arcsin(tau * np.sin(t))
    x = 2 * k * np.sin(th_cm / 2) * R
    F = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    print(f"theta_lab = {th_lab:5.1f}°  theta_cm = {np.degrees(th_cm):5.1f}°  qR = {x:4.2f}  |F|^2 = {F**2:.3f}  medido = {r:.2f}")
