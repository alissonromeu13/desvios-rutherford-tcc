"""Figura 2 do TCC (Seção 2.1.2) – Anel de incidência e cone de saída. Partículas com parâmetro de impacto entre
b e b + db (anel de área dσ = 2πb db) são espalhadas no cone entre θ e θ − |dθ| (ângulo sólido dΩ = 2π sen θ dθ).
Trajetórias calculadas por u(φ) = (D/2b²)(cos φ − 1) + sen φ / b (Seção 2.1.1). Figura original do autor."""
import numpy as np
from estilo import *
from matplotlib.patches import Ellipse, Arc

D = 0.8
th1, th2 = np.radians(62.0), np.radians(46.0)          # b -> θ1 ; b + db -> θ2 < θ1
b1, b2 = D / 2 / np.tan(th1 / 2), D / 2 / np.tan(th2 / 2)
def traj(b, th):
    ph = np.linspace(0.01, np.pi - th - 0.01, 3000)
    u = D / (2 * b**2) * (np.cos(ph) - 1) + np.sin(ph) / b
    r = 1 / u; return -r * np.cos(ph), r * np.sin(ph)
fig, ax = plt.subplots(figsize=(6.4, 3.9)); ax.set_aspect("equal"); ax.axis("off")
Rs = 4.6
ax.plot([-5.2, Rs + 0.4], [0, 0], color=TINTA, lw=0.7, ls=(0, (6, 3)))
ax.plot(0, 0, "o", color=TINTA, ms=6); ax.text(0.1, -0.35, r"$+Ze$", fontsize=9.5)
for b, th, c in [(b1, th1, AZUL), (b2, th2, LARANJA)]:
    X, Y = traj(b, th); m = (X > -5.0) & (np.hypot(X, Y) < Rs)
    ax.plot(X[m], Y[m], color=c, lw=1.5)
    ax.plot(X[m][::-1][:1], Y[m][::-1][:1], "o", color=c, ms=2.5)
# anel de entrada (perspectiva) em x = -4.3
xa = -4.3
for b, c in [(b1, AZUL), (b2, LARANJA)]:
    ax.add_patch(Ellipse((xa, 0), 0.6 * b, 2 * b, fill=False, color=c, lw=1.0))
ax.annotate("", xy=(xa - 0.5, b1), xytext=(xa - 0.5, 0), arrowprops=dict(arrowstyle="|-|,widthA=0.25,widthB=0.25", color=TINTA, lw=0.8, shrinkA=0, shrinkB=0))
ax.text(xa - 0.58, b1 / 2, r"$b$", ha="right", va="center", fontsize=10)
ax.annotate("", xy=(xa - 0.95, b2), xytext=(xa - 0.95, b1), arrowprops=dict(arrowstyle="|-|,widthA=0.25,widthB=0.25", color=TINTA, lw=0.8, shrinkA=0, shrinkB=0))
ax.text(xa - 1.03, (b1 + b2) / 2, r"$\mathrm{d}b$", ha="right", va="center", fontsize=9)
ax.text(xa, -b2 - 0.25, r"$\mathrm{d}\sigma = 2\pi b\,\mathrm{d}b$", ha="center", va="top", fontsize=9.5)
# esfera de raio Rs e faixa (anel) entre θ2 e θ1 em perspectiva
ax.add_patch(Arc((0, 0), 2 * Rs, 2 * Rs, theta1=5, theta2=120, color=CINZA, lw=0.7, ls=":"))
for th, c in [(th1, AZUL), (th2, LARANJA)]:
    ax.plot([0, Rs * np.cos(th)], [0, Rs * np.sin(th)], color=c, lw=0.6, ls=(0, (2, 2)))
    xc = Rs * np.cos(th); h = Rs * np.sin(th)
    ax.add_patch(Arc((xc, 0), 0.5 * h, 2 * h, theta1=90, theta2=270, color=c, lw=0.7, alpha=0.8))
    ax.add_patch(Arc((xc, 0), 0.5 * h, 2 * h, theta1=-90, theta2=90, color=c, lw=0.6, ls=(0, (2, 2))))
tt = np.linspace(th2, th1, 50); ax.plot(Rs * np.cos(tt), Rs * np.sin(tt), color=TINTA, lw=2.2)
ax.text(Rs * np.cos((th1 + th2) / 2) + 0.12, Rs * np.sin((th1 + th2) / 2) + 0.12, r"$\mathrm{d}\Omega = 2\pi\,\mathrm{sen}\,\theta\,\mathrm{d}\theta$", fontsize=9.5)
ax.add_patch(Arc((0, 0), 1.6, 1.6, theta1=0, theta2=np.degrees(th1), color=AZUL, lw=0.8))
ax.text(0.95 * np.cos(th1 / 2), 0.95 * np.sin(th1 / 2) - 0.08, r"$\theta$", fontsize=10, color=AZUL)
ax.add_patch(Arc((0, 0), 2.5, 2.5, theta1=np.degrees(th2), theta2=np.degrees(th1), color=TINTA, lw=0.8))
ax.text(1.42 * np.cos((th1 + th2) / 2), 1.42 * np.sin((th1 + th2) / 2), r"$|\mathrm{d}\theta|$", fontsize=9)
ax.text(Rs + 0.1, 0.12, r"$r$", fontsize=9, color=CINZA)
ax.set(xlim=(-5.4, 5.9), ylim=(-Rs * 0.95, Rs * 1.02))
salvar(fig, "2/fig_rutherford_feixe")
