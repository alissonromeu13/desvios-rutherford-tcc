"""Fator de forma de carga confrontado com Bieler (Seção 4.1.2).
Na primeira aproximação de Born, sigma/sigma_R = |F(q)|^2. Esfera homogênea, R = 1,2 A^(1/3) fm,
partículas do RaC no alumínio com E_lab = 7,35 MeV (Apêndice D), k = 1,033 fm^-1 no centro de massa,
ângulos efetivos de laboratório dos quatro diafragmas (Tabela do Apêndice D) convertidos ao centro de massa (Apêndice C).
Também as distribuições exponencial e gaussiana, com o mesmo raio quadrático médio: a conclusão não depende do modelo."""
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

# Mesma comparação com as distribuições exponencial e gaussiana da Seção 2.4, com o mesmo raio quadrático médio
# da esfera, <r^2> = 3R^2/5: exponencial, F = (1 + q^2 a^2)^-2 com <r^2> = 12 a^2; gaussiana, F = exp(-q^2 <r^2>/6).
r2 = 0.6 * R**2
for th_lab, r in zip((20.2, 30.9, 50.7, 75.8), medido):
    t = np.radians(th_lab)
    q = 2 * k * np.sin((t + np.arcsin(tau * np.sin(t))) / 2)
    Fx = (1 + q**2 * r2 / 12) ** -2
    Fg = np.exp(-q**2 * r2 / 6)
    print(f"theta_lab = {th_lab:5.1f}°  |F|^2 exponencial = {Fx**2:.3f}  gaussiana = {Fg**2:.3f}  medido = {r:.2f}")

# Modelo da Seção 4.1.2, sigma/sigma_R = [1 - g x j1(x)]^2: x j1(x) tem máximo em x = 2,744;
# o mínimo da razão fica nesse ponto só se g x j1(x) < 1 ali, isto é, g < 1/max.
from scipy.optimize import minimize_scalar
xj1 = lambda x: np.sin(x) / x - np.cos(x)
r = minimize_scalar(lambda x: -xj1(x), bounds=(1, 4), method="bounded")
print(f"max de x j1(x) = {-r.fun:.4f} em x = {r.x:.3f}; o mínimo em qR = 2,74 exige g < {1 / -r.fun:.2f} (valores do texto: 0,35 e 0,21)")
