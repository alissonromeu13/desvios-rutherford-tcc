"""Painel (b) retirado da Figura 13 da versão 1.1 (disco negro com borda suave); guardado em extras/. Origem da escolha Delta = 0,5.
S_l = 1/[1 + exp((kR - l - 1/2)/Delta)];  Delta = 0 é o corte abrupto (S_l = 0 para l + 1/2 < kR).
Com o corte abrupto, sigma_reação/(pi R^2) = sum_l (2l+1)(1 - |S_l|^2)/(kR)^2 salta cada vez que kR + 1/2 cruza
um inteiro (o "dente de serra" do l discreto). O script mede a oscilação em torno da tendência suave para kR > 9
e o valor em kR = 10 (o desvio de ordem Delta/kR discutido no texto da Seção 4.2.4.2).
Resultado: corte abrupto ~ 0,106 (cerca de 10 % de pi R^2); Delta = 0,3 -> 0,005; Delta = 0,5 -> 0,0003 (< 0,1 %)."""
import numpy as np

kR = np.linspace(5, 40, 3501)

def sigma_r(D):
    out = []
    for x in kR:
        l = np.arange(int(x + 15 * max(D, 0.01) + 5))
        S = (l + 0.5 > x).astype(float) if D == 0 else 1 / (1 + np.exp((x - l - 0.5) / D))
        out.append(np.sum((2 * l + 1) * (1 - S**2)) / x**2)
    return np.array(out)

for D in (0.0, 0.1, 0.2, 0.3, 0.5, 0.7):
    s = sigma_r(D)
    tendencia = np.polyval(np.polyfit(1 / kR, s, 3), 1 / kR)
    m = kR > 9
    print(f"Delta = {D:.1f}: oscilação máxima (kR > 9) = {np.max(np.abs(s - tendencia)[m]):.4f} pi R^2;"
          f"  sigma_r(kR = 10) = {s[np.argmin(abs(kR - 10))]:.3f} pi R^2")
