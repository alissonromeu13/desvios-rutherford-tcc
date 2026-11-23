"""Raio de corte da condição de fase (Seção 2.3.3).
Ao longo de uma reta com parâmetro de impacto b, a fase acumulada no potencial coulombiano blindado
V = k_e zZe^2 e^(-r/R)/r é (1/hbar v) int |V| ds = 2 eta K0(b/R)  (Abramowitz e Stegun, p. 376, eq. 9.6.24).
A condição 2 eta K0(r_lim/R) = 1 é resolvida numericamente e comparada às formas assintóticas:
  eta << 1:  r_lim ~ 2R exp(-gamma_E) exp(-1/(2 eta))       (eq. 9.6.13, p. 375)
  eta >> 1:  x = r_lim/R ~ ln(sqrt(2 pi) eta) - (1/2) ln ln(sqrt(2 pi) eta)   (eq. 9.7.2, p. 378)
Também imprime o critério de altas energias de Sakurai (p. 387, eq. 6.85), 2 eta ln(k a) << 1, e a fase na escala
do comprimento de onda reduzido para os elétrons de 10 keV no hidrogênio."""
import numpy as np
from scipy.special import k0
from scipy.optimize import brentq

gE = 0.5772156649
alfa = 1 / 137.035999084
fase = lambda x, eta: 2 * eta * k0(x)
raiz = lambda eta: brentq(lambda x: fase(x, eta) - 1, 1e-14, 60)

print("alfa + H (eta = 0,23):  numérico r_lim/R = %.3f;  fórmula 2 e^-gamma e^(-1/2eta) = %.3f"
      % (raiz(0.23), 2 * np.exp(-gE) * np.exp(-1 / 0.46)))
eta = 18.36   # ouro a 7,35 MeV, a energia das partículas de Bieler (Seção 4.1.1)
y = np.log(np.sqrt(2 * np.pi) * eta)
print("ouro (eta = 18,4):      numérico r_lim/R = %.2f;  assintótica = %.2f;  ln(2 eta) = %.2f" % (raiz(eta), y - 0.5 * np.log(y), np.log(2 * eta)))

# elétrons de 10 keV no hidrogênio
me = 0.51099895; T = 0.010
b = np.sqrt(1 - 1 / (1 + T / me) ** 2); eta_e = alfa / b
lam = 197.3269804 / (np.sqrt(T * (T + 2 * me))) * 1e-5      # hbar/p em Å (1 fm = 1e-5 Å)
rB = 0.529177
print("elétrons de 10 keV: beta = %.3f, eta = %.4f, lambda_bar = %.4f Å = %.4f r_B, k r_B = %.1f"
      % (b, eta_e, lam, lam / rB, rB / lam))
print("   critério de altas energias: 2 eta ln(k r_B) = %.3f" % (2 * eta_e * np.log(rB / lam)))
print("   r_lim/r_B: numérico = %.2e;  fórmula = %.2e" % (raiz(eta_e), 2 * np.exp(-gE) * np.exp(-1 / (2 * eta_e))))
print("   fase em b = lambda_bar: 2 eta K0(lambda_bar/r_B) = %.3f rad" % fase(lam / rB, eta_e))
