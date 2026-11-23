"""Figura 1 do TCC (Seção 2.1.1). Vetores em negrito, como no texto – Trajetória hiperbólica de Rutherford. Figura original do autor.
Notação do texto: u = 1/r = (D/2b²)(cos φ − 1) + sen φ / b, com φ medido a partir do eixo que aponta para a
partícula incidente (−x); entrada em φ = 0, saída em φ = π − θ, cot(θ/2) = 2b/D, b = b',
r_min = (D/2)[1 + 1/sen(θ/2)] em φ = (π − θ)/2.
No ponto do ramo de saída: φ' é o ângulo polar e β = π − φ' − θ é o ângulo entre r' e v' (b' = r' sen β)."""
import numpy as np
from estilo import *
from matplotlib.patches import Arc, FancyArrowPatch

D = 1.0
th = np.radians(65.0)
b = D / 2 / np.tan(th / 2)
phi_max = np.pi - th
def r_of(phi):
    u = D / (2 * b**2) * (np.cos(phi) - 1) + np.sin(phi) / b
    return 1 / u
def xy(phi):
    r = r_of(phi); return -r * np.cos(phi), r * np.sin(phi)
ph = np.linspace(0.012, phi_max - 0.012, 3000)
X, Y = xy(ph)
m = (X > -4.2) & (X**2 + Y**2 < 4.6**2)
rmin = D / 2 * (1 + 1 / np.sin(th / 2)); ph_min = (np.pi - th) / 2
print(f"b = {b:.3f} D, r_min = {rmin:.3f} D")

fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.set_aspect("equal"); ax.axis("off")
ax.plot([-4.3, 2.6], [0, 0], color=TINTA, lw=0.8)                         # eixo
ax.plot([-4.3, 1.2], [b, b], color=CINZA, lw=0.8, ls="--")                 # assíntota de entrada
d = np.array([np.cos(th), np.sin(th)]); n = np.array([np.sin(th), -np.cos(th)])  # direção e normal da saída
p0 = -b * n                                                                  # assíntota de saída passa a b' de O
s = np.linspace(0.6, 4.6, 2); ax.plot(p0[0] + s * d[0], p0[1] + s * d[1], color=CINZA, lw=0.8, ls="--")
ax.plot([0, 3.7 * d[0]], [0, 3.7 * d[1]], color=CINZA, lw=0.7, ls=":")      # paralela por O (define θ)
ax.plot(X[m], Y[m], color=AZUL, lw=1.8)
i = np.argmin(abs(X[m] + 2.6)); ax.annotate("", xy=(X[m][i+40], Y[m][i+40]), xytext=(X[m][i], Y[m][i]),
    arrowprops=dict(arrowstyle="-|>", color=AZUL, lw=1.2))
j = np.argmin(abs(np.hypot(X[m], Y[m]) - 3.6) + (Y[m] < 1) * 9); ax.annotate("", xy=(X[m][j+40], Y[m][j+40]), xytext=(X[m][j], Y[m][j]),
    arrowprops=dict(arrowstyle="-|>", color=AZUL, lw=1.2))
ax.plot(0, 0, "o", color=TINTA, ms=7); ax.text(0.12, -0.28, r"$+Ze$", fontsize=10)
# b e b'
xb = -3.7; ax.annotate("", xy=(xb, b), xytext=(xb, 0), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=0.8))
ax.text(xb + 0.08, b / 2, r"$b$", fontsize=10, va="center")
ax.annotate("", xy=(p0[0] + 3.4*d[0], p0[1] + 3.4*d[1]), xytext=(3.4*d[0], 3.4*d[1]), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=0.8))
mid = 3.4 * d + p0 / 2; ax.text(mid[0] - 0.05, mid[1] + 0.12, r"$b'$", fontsize=10)
# θ
ax.add_patch(Arc((0, 0), 1.5, 1.5, theta1=0, theta2=np.degrees(th), color=TINTA, lw=0.8))
ax.text(0.85 * np.cos(th / 2), 0.85 * np.sin(th / 2), r"$\theta$", fontsize=11, ha="left", va="center")
# ponto genérico P: r e φ
phP = 0.36; xP, yP = xy(phP)
ax.annotate("", xy=(xP, yP), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=LARANJA, lw=1.1))
ax.text(xP * 0.45, yP * 0.45 - 0.2, r"$\mathbf{r}$", fontsize=10, color=LARANJA, ha="right")
ax.add_patch(Arc((0, 0), 2.6, 2.6, theta1=180 - np.degrees(phP), theta2=180, color=LARANJA, lw=0.8))
ax.text(-1.42 * np.cos(phP / 2), 1.42 * np.sin(phP / 2) - 0.02, r"$\varphi$", fontsize=10, color=LARANJA, ha="right", va="center")
# velocidade em P e componentes
er = np.array([xP, yP]) / np.hypot(xP, yP); ep = np.array([-er[1], er[0]]) * -1
dphi = 1e-4; x2, y2 = xy(phP + dphi); t = np.array([x2 - xP, y2 - yP]); t /= np.hypot(*t)
L = 0.95; ax.annotate("", xy=(xP + L * t[0], yP + L * t[1]), xytext=(xP, yP), arrowprops=dict(arrowstyle="-|>", color=TINTA, lw=1.0))
ax.text(xP + L * t[0], yP + L * t[1] + 0.08, r"$\mathbf{v}$", fontsize=10, ha="center", va="bottom")
vr = np.dot(t, er) * L * er; vp = t * L - vr
ax.annotate("", xy=(xP + vr[0], yP + vr[1]), xytext=(xP, yP), arrowprops=dict(arrowstyle="-|>", color=CINZA, lw=0.8))
ax.text(xP + vr[0] - 0.05, yP + vr[1] + 0.06, r"$v_r$", fontsize=8.5, color=CINZA)
ax.annotate("", xy=(xP + vp[0], yP + vp[1]), xytext=(xP, yP), arrowprops=dict(arrowstyle="-|>", color=CINZA, lw=0.8))
ax.text(xP + vp[0] - 0.08, yP + vp[1] + 0.02, r"$v_\varphi$", fontsize=8.5, color=CINZA, ha="right")
ax.plot(xP, yP, "o", color=TINTA, ms=3.5)
# ponto P' no ramo de saída: r', φ', v' e β
phQ = phi_max - 0.55; xQ, yQ = xy(phQ)
ax.annotate("", xy=(xQ, yQ), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=LARANJA, lw=1.1))
ax.text(xQ * 0.62 + 0.12, yQ * 0.62, r"$\mathbf{r}'$", fontsize=10, color=LARANJA, ha="left", va="center")
ax.add_patch(Arc((0, 0), 0.9, 0.9, theta1=180 - np.degrees(phQ), theta2=180, color=LARANJA, lw=0.8))
ax.text(0.66 * np.cos(np.radians(140)), 0.66 * np.sin(np.radians(140)), r"$\varphi'$", fontsize=9, color=LARANJA, ha="center", va="center")
x2, y2 = xy(phQ + 1e-4); tq = np.array([x2 - xQ, y2 - yQ]); tq /= np.hypot(*tq)
ax.annotate("", xy=(xQ + L * tq[0], yQ + L * tq[1]), xytext=(xQ, yQ), arrowprops=dict(arrowstyle="-|>", color=TINTA, lw=1.0))
ax.text(xQ + L * tq[0] + 0.08, yQ + L * tq[1], r"$\mathbf{v}'$", fontsize=10, ha="left", va="center")
ax.plot([xQ, xQ + 0.75 * xQ / np.hypot(xQ, yQ)], [yQ, yQ + 0.75 * yQ / np.hypot(xQ, yQ)], color=LARANJA, lw=0.7, ls=":")
aR = np.degrees(np.arctan2(yQ, xQ)); aV = np.degrees(np.arctan2(tq[1], tq[0]))
ax.add_patch(Arc((xQ, yQ), 1.1, 1.1, theta1=min(aR, aV), theta2=max(aR, aV), color=TINTA, lw=0.8))
am = np.radians((aR + aV) / 2); ax.text(xQ + 0.68 * np.cos(am), yQ + 0.68 * np.sin(am), r"$\beta$", fontsize=10, ha="center", va="center")
ax.plot(xQ, yQ, "o", color=TINTA, ms=3.5)
# máxima aproximação e eixo de simetria
xm, ym = xy(ph_min)
ax.plot([0, 1.45 * xm], [0, 1.45 * ym], color=CINZA, lw=0.7, ls="-.")
ax.annotate("", xy=(xm, ym), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=VERDE, lw=1.1))
ax.text(0.8 * xm - 0.08, 0.8 * ym, r"$r_{\text{min}}$", fontsize=10, color=VERDE, ha="right", va="center")
ax.plot(xm, ym, "o", color=TINTA, ms=3.5)
ax.text(1.45 * xm - 0.05, 1.45 * ym + 0.05, "eixo de\nsimetria", fontsize=7, color=CINZA, ha="right", va="bottom")
ax.text(-4.25, b + 0.12, r"$\varphi \to 0$", fontsize=8.5, color=CINZA)
ax.text(-4.25, b + 0.45, "partícula alfa", fontsize=8.5, color=TINTA, ha="left", va="bottom")
e = p0 + 4.1 * d; ax.text(e[0] - 0.35, e[1], r"$\varphi \to \pi - \theta$", fontsize=8.5, color=CINZA, ha="right")
ax.set(xlim=(-4.4, 2.9), ylim=(-0.45, 4.3))
salvar(fig, "2/fig_rutherford_trajetoria")
