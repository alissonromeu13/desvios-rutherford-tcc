"""Origem da escolha K0R = 5,30 na Figura 12 e no diagrama de Argand de extras/ (poço quadrado atrativo, unidades de R).
Varre K0R e registra, para l = 1, 2 e 3, a posição (kR) e a largura aproximada (faixa de kR com sin^2 delta_l > 1/2)
do primeiro máximo de sin^2 delta_l com kR < 2,5. Critério usado: uma ressonância de forma estreita, em energia
baixa o bastante para que poucas ondas parciais contribuam, e larga o bastante para ser resolvida no gráfico.
Com K0R = 5,3, a onda l = 3 ressoa em kR ≈ 1,83 com largura ≈ 0,12 (para K0R = 5,6 a largura já cai a 0,01).
Fórmula de tan delta_l: a mesma da Figura 12 (Seção 4.2.3)."""
import numpy as np
from scipy.special import spherical_jn as j, spherical_yn as y

def tand(l, kR, K0R):
    KR = np.sqrt(kR**2 + K0R**2)
    num = kR * j(l, kR, derivative=True) * j(l, KR) - KR * j(l, kR) * j(l, KR, derivative=True)
    den = kR * y(l, kR, derivative=True) * j(l, KR) - KR * y(l, kR) * j(l, KR, derivative=True)
    return num / den

kR = np.linspace(0.01, 3, 3000)
for K0R in np.round(np.arange(3.0, 6.0, 0.1), 1):
    for l in (1, 2, 3):
        s = np.sin(np.arctan(tand(l, kR, K0R))) ** 2
        i = np.argmax(s)
        if s[i] > 0.99 and kR[i] < 2.5:
            print(f"K0R = {K0R:.1f}  l = {l}: ressonância em kR = {kR[i]:.2f}, largura ≈ {np.sum(s > 0.5) * (kR[1] - kR[0]):.3f}")
