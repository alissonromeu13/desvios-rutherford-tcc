"""Figura 14 do TCC (Seção 4.2.4.2) – Formação e desaparecimento da sombra do disco negro, e o padrão de Airy.

Onda plana incidente e^{ikz} sobre um disco opaco de raio R. Pelo princípio de Babinet, o campo atrás do
disco é o incidente menos o campo transmitido por uma abertura de mesmo raio:  psi = 1 - psi_ab.
Aproximação de Fresnel (paraxial; válida para kR >> 1 e z >> R):
    psi_ab(rho, z) = -i u e^{i u rho~^2/2} int_0^1 J0(u rho~ s) e^{i u s^2/2} s ds,
    u = k R^2 / z = 2 pi N_F,   N_F = R^2/(lambda z) (número de Fresnel),   rho~ = rho/R.
As imagens dependem só de N_F. No eixo, |psi|^2 = 1 exatamente (mancha de Poisson-Arago).
Campo distante da onda espalhada (Fraunhofer): intensidade relativa [2 J1(x)/x]^2, x = k R theta,
primeiro zero em theta_1 = 3,8317/(kR)  ->  mesmo padrão de Airy da Eq. de d sigma/d Omega da seção.
"""
import numpy as np
from scipy.special import j0, j1, jn_zeros
from estilo import *
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.ticker import FixedLocator, FixedFormatter
import os

def psi_disco(rho, u, M=8000):
    s = (np.arange(M) + 0.5) / M
    integ = (j0(u * np.outer(rho, s)) * np.exp(1j * u * s**2 / 2) * s).sum(1) / M
    return 1 - (-1j * u * np.exp(1j * u * rho**2 / 2) * integ)

# ------------------------------------------------------------------ checagens
for NF in (10, 1, 0.01):
    I0 = abs(psi_disco(np.array([0.0]), 2 * np.pi * NF))**2
    print(f"   N_F = {NF}: |psi(0)|^2 = {I0[0]:.6f}"); assert abs(I0 - 1) < 1e-3
print("mancha de Poisson-Arago: |psi(0)|^2 = 1 em todas as distâncias")
u = 0.05; r = np.linspace(0.5, 120, 7)
razao = np.abs(1 - psi_disco(r, u, 20000)) / np.abs(u * j1(u * r) / (u * r))
assert np.allclose(razao, 1, atol=0.02)
print("limite de Fraunhofer: onda espalhada coincide com o padrão de Airy (desvio < 2%)")
x1, x2 = jn_zeros(1, 2)
lam_fm, R_fm = 6.42, 7.0
print(f"exemplo nuclear (alfa de 5 MeV em Au): a sombra só existe até z ~ R^2/lambda = {R_fm**2/lam_fm:.1f} fm; "
      f"num detector a 10 cm, N_F ~ {R_fm**2/(lam_fm*1e14):.0e}")

# ------------------------------------------------------------------ figura
W, H = 7.2, 4.0
fig = plt.figure(figsize=(W, H))
box = lambda x0, y0, w, h: fig.add_axes([x0 / W, y0 / H, w / W, h / H])
S = 1.25                                                       # lado das imagens (polegadas)

# (topo) esquema fora de escala
ax0 = box(0.05, 1.98, 7.1, 1.98)
ax0.set_xlim(-3.2, 15.2); ax0.set_ylim(-2.6, 2.9); ax0.axis("off")
for zw in np.arange(-2.8, -0.8, 0.45):
    ax0.plot([zw, zw], [-1.9, 1.9], color=AZUL, lw=0.8)
ax0.add_patch(FancyArrowPatch((-2.9, 2.35), (-0.9, 2.35), arrowstyle="->", mutation_scale=9, color=AZUL, lw=0.9))
ax0.text(-1.85, 2.5, "onda plana incidente", fontsize=7, color=AZUL, ha="center", va="bottom")
grad = np.zeros((2, 256, 4)); grad[..., :3] = 0.42
grad[..., 3] = np.linspace(0.45, 0.0, 256)[None, :]             # sombra que se desfaz com a distância
ax0.imshow(grad, extent=(0, 7.5, -1, 1), aspect="auto", interpolation="bilinear", zorder=1)
ax0.text(2.4, 0, "sombra", fontsize=7, color=TINTA, ha="center", va="center", zorder=3)
th_s = 0.17                                                    # semiabertura do cone (esquemática)
zc = np.array([0, 13.6])
for sgn in (1, -1):
    ax0.plot(zc, sgn * np.tan(th_s) * zc, color=LARANJA, lw=0.9, ls="--", zorder=2)
ax0.add_patch(Rectangle((-0.09, -1), 0.18, 2, color="black", zorder=4))
ax0.annotate("", xy=(0.45, 1.0), xytext=(0.45, 0.0), arrowprops=dict(arrowstyle="<->", lw=0.6, mutation_scale=6))
ax0.text(0.55, 0.5, "$R$", fontsize=7.5, va="center")
ax0.text(0, -1.25, "disco negro", fontsize=7, ha="center", va="top")
for z, lab in [(1.4, "(a)"), (5.6, "(b)"), (9.8, "(c)")]:
    ax0.plot([z, z], [-2.0, 2.0], color=TINTA, lw=1.0, zorder=3)
    ax0.text(z, 2.1, r"\textbf{" + lab + "}", fontsize=7.5, ha="center", va="bottom")
ax0.plot([13.6, 13.6], [-2.4, 2.4], color=TINTA, lw=1.4)
ax0.text(13.6, 2.5, r"\textbf{(d)} anteparo distante", fontsize=7, ha="center", va="bottom")
ax0.text(11.7, 0.0, "cone de difração\n" + r"$\theta_1 \approx 3{,}83/kR$", fontsize=7, color=LARANJA, ha="center", va="center")
ax0.add_patch(FancyArrowPatch((0.3, -2.45), (13.3, -2.45), arrowstyle="->", mutation_scale=8, color=CINZA, lw=0.8))
ax0.text(6.8, -2.38, r"distância $z$ crescente: $N_F = R^2/\lambda z$ diminui (esquema fora de escala)",
         fontsize=6.8, color=CINZA, ha="center", va="bottom")

def moldura(ax):
    ax.grid(False)
    for s_ in ax.spines.values():
        s_.set_visible(True); s_.set_linewidth(0.5)
    ax.tick_params(labelsize=7, length=2)

# (a)-(c) intensidade |psi|^2 em anteparos a distâncias crescentes
L, N = 2.5, 501
xg = np.linspace(-L, L, N); X, Y = np.meshgrid(xg, xg); Rg = np.hypot(X, Y)
rr = np.linspace(0, L * np.sqrt(2) + 0.01, 1400)
vmax = 1.3
for i, (NF, lab) in enumerate([(10, "a"), (1, "b"), (0.01, "c")]):
    ax = box(0.42 + i * (S + 0.10), 0.42, S, S)
    Ir = np.abs(psi_disco(rr, 2 * np.pi * NF)) ** 2
    I2 = np.interp(Rg, rr, Ir)
    im = ax.imshow(I2, extent=(-L, L, -L, L), origin="lower", cmap="gray", vmin=0, vmax=vmax, interpolation="bilinear")
    ax.add_patch(Circle((0, 0), 1, fill=False, ec=LARANJA, lw=0.7, ls="--"))
    ax.set_xticks([-2, 0, 2]); ax.set_yticks([-2, 0, 2])
    ax.set_xticklabels(["$-2$", "0", "2"]); ax.set_yticklabels(["$-2$", "0", "2"] if i == 0 else [])
    moldura(ax)
    ax.set_xlabel("$x/R$", fontsize=8, labelpad=1)
    if i == 0:
        ax.set_ylabel("$y/R$", fontsize=8, labelpad=1)
    nf = f"{NF:g}".replace(".", "{,}")
    ax.set_title(r"\textbf{(" + lab + r")}" + f" $N_F = {nf}$", fontsize=8, pad=3, loc="left")
    print(f"N_F = {NF}: I(centro) = {Ir[0]:.3f}, I mín = {Ir.min():.3f}, I(rho = 2R) = {np.interp(2, rr, Ir):.3f}")
cax = box(0.42 + 3 * S + 0.27, 0.42, 0.08, S)
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([0, 0.5, 1]); cb.set_ticklabels(["0", "0,5", "1"])
cb.ax.tick_params(labelsize=6.5, length=2)
cb.set_label(r"$|\psi|^2/|\psi_\mathrm{inc}|^2$", fontsize=7.5, labelpad=2)

# (d) padrão de Airy da onda espalhada (campo distante)
axd = box(5.35, 0.42, S, S)
T = 3.5
tg = np.linspace(-T, T, 601); TX, TY = np.meshgrid(tg, tg)
xx = np.hypot(TX, TY) * x1
xx[xx == 0] = 1e-12
Iff = (2 * j1(xx) / xx) ** 2
imd = axd.imshow(np.log10(Iff + 1e-12), extent=(-T, T, -T, T), origin="lower", cmap="gray", vmin=-4, vmax=0,
                 interpolation="bilinear")
for rad in (1, x2 / x1):
    axd.add_patch(Circle((0, 0), rad, fill=False, ec=LARANJA, lw=0.6, ls="--"))
axd.text(1.12, 0.05, r"$\theta_1$", fontsize=7.5, color=LARANJA, va="center")
axd.set_xticks([-3, 0, 3]); axd.set_yticks([-3, 0, 3])
axd.set_xticklabels(["$-3$", "0", "3"]); axd.set_yticklabels(["$-3$", "0", "3"])
moldura(axd)
axd.set_xlabel(r"$\theta_x/\theta_1$", fontsize=8, labelpad=1)
axd.set_ylabel(r"$\theta_y/\theta_1$", fontsize=8, labelpad=0)
axd.set_title(r"\textbf{(d)} onda espalhada", fontsize=8, pad=3, loc="left")
caxd = box(5.35 + S + 0.07, 0.42, 0.08, S)
cbd = fig.colorbar(imd, cax=caxd)
cbd.set_ticks([-4, -2, 0]); cbd.set_ticklabels(["$10^{-4}$", "$10^{-2}$", "1"])
cbd.ax.tick_params(labelsize=6.5, length=2)
cbd.set_label("intensidade relativa", fontsize=7.5, labelpad=2)
salvar(fig, "6.3/fig12_sombra_disco_negro", bbox_inches=None)
