"""Figura 8 do TCC (Seção 2.4.4) – Validade da 1ª aproximação de Born: potencial gaussiano V(r) = V0 exp(-r^2/a^2)
(o mesmo perfil da Seção 2.4.3). Exato = ondas parciais (Numerov); Born = forma fechada.

  f_B(q)     = -(mu V0 sqrt(pi) a^3 / (2 hbar^2)) exp(-q^2 a^2/4) = (xi sqrt(pi) a / 4) exp(-q^2 a^2/4)   (V0 < 0)
  sigma_B(E) = (pi^2/8) xi^2 a^2 [1 - exp(-2 k^2 a^2)] / (k^2 a^2),      sigma_B(0) = pi^2 xi^2 a^2 / 4
  xi         = 2 mu |V0| a^2 / hbar^2   (Born exige xi << 1 em baixas energias)
Comprimento de espalhamento: exato a0 = -lim tan(delta0)/k ; Born a_B = -f_B(0) = -xi sqrt(pi) a / 4.
"""
import numpy as np
from estilo import *
from scatter_core import smatrix, nuclear_amplitude, HBARC, M_N

mu, a = M_N, 2.0
casos = [(-2.0, AZUL), (-20.0, LARANJA)]
xi_de = lambda V0: 2 * mu * abs(V0) * a**2 / HBARC**2
fb = lambda V0, q: (xi_de(V0) * np.sqrt(np.pi) * a / 4) * np.exp(-(q * a) ** 2 / 4) * np.sign(-V0)

def sig_born(V0, E):
    """Seção de choque total de Born em forma fechada (fm^2)."""
    y = 2 * mu * E / HBARC**2 * a**2          # k^2 a^2
    return np.pi**2 / 8 * xi_de(V0) ** 2 * a**2 * (-np.expm1(-2 * y)) / y

def exato(V0, E, h=0.005):
    k = np.sqrt(2 * mu * E) / HBARC
    lmax = int(k * 3 * a) + 8
    return smatrix(lambda r: V0 * np.exp(-(r / a) ** 2) + 0j, mu, E, lmax, 6 * a, h=h)

def nos_E0(V0, h=0.002):
    """Número de nós da solução de onda s em E = 0 (= número de estados ligados s)."""
    r = np.arange(1, int(8 * a / h)) * h
    g = -2 * mu * V0 * np.exp(-(r / a) ** 2) / HBARC**2
    w = 1 + h * h * g / 12; u = np.zeros_like(r); u[0], u[1] = r[0], r[1]
    for i in range(1, len(r) - 1):
        u[i + 1] = ((12 - 10 * w[i]) * u[i] - w[i - 1] * u[i - 1]) / w[i + 1]
    return int(np.sum(np.diff(np.sign(u)) != 0))

# ---- checagem da forma fechada contra integração numérica
for V0, _ in casos:
    for E in (0.01, 20.0, 300.0):
        k = np.sqrt(2 * mu * E) / HBARC
        x = np.linspace(-1, 1, 200001); q = k * np.sqrt(2 - 2 * x)
        num = 2 * np.pi * np.trapezoid(fb(V0, q) ** 2, x)
        assert abs(num / sig_born(V0, E) - 1) < 1e-6, (V0, E)
print("forma fechada de sigma_B confere com a integral numérica (erro < 1e-6)")

E = np.logspace(-2, np.log10(300), 50)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.6))
handles = []
info = {}
for V0, c in casos:
    xi = xi_de(V0)
    sx, opt = [], []
    for e in E:
        S, k = exato(V0, e)
        l = np.arange(len(S))
        st = np.pi / k**2 * np.sum((2 * l + 1) * np.abs(1 - S) ** 2)
        f0 = nuclear_amplitude(S, 0.0, k, np.array([0.0]))[0]
        opt.append(abs(4 * np.pi / k * f0.imag - st) / st)
        sx.append(st * 10)
    sb = sig_born(V0, E) * 10
    S, k = exato(V0, 1e-4, h=0.002)
    a0 = -np.tan(np.angle(S[0]) / 2) / k
    aB = -fb(V0, 0.0)
    nb = nos_E0(V0)
    info[V0] = (xi, a0, aB, nb)
    print(f"V0={V0}: xi={xi:.3f}; estados ligados s = {nb}; a0 exato = {a0:+.3f} fm, a_B = {aB:+.3f} fm; "
          f"sigma(0): exato {4*np.pi*a0**2*10:.0f} mb, Born {np.pi**2*xi**2*a**2/4*10:.0f} mb "
          f"(razão {a0**2/aB**2:.2f}); teorema óptico erro máx = {max(opt):.1e}; Born/exato em 300 MeV = {sb[-1]/sx[-1]:.3f}")
    rot = f"$V_0 = " + f"{V0:.1f}".replace(".", "{,}") + "$ MeV ($\\xi = " + f"{xi:.2f}".replace(".", "{,}") + "$)"
    h1, = a1.plot(E, sx, color=c, ls="-", label=rot + " – exato")
    h2, = a1.plot(E, sb, color=c, ls="--", label=rot + " – Born")
    handles += [h1, h2]

xi, a0, aB, nb = info[-20.0]
txt = ("estado ligado raso:\n"
       + f"$a_0 = {a0:+.1f}$ fm (exato)\n".replace(".", "{,}")
       + f"$a_B = {aB:+.1f}$ fm (Born)".replace(".", "{,}"))
a1.annotate(txt, xy=(0.012, 4 * np.pi * a0**2 * 10 * 0.97), xytext=(0.03, 25), fontsize=7, color=LARANJA,
            arrowprops=dict(arrowstyle="->", color=LARANJA, lw=0.7))
a1.set(xscale="log", yscale="log", xlabel=r"$E_\mathrm{cm}$ [MeV]", ylabel=r"$\sigma$ [mb]", ylim=(0.05, 2e4))

Ed = 20.0
th = np.linspace(0.001, np.pi, 400)
for V0, c in casos:
    S, k = exato(V0, Ed)
    f = nuclear_amplitude(S, 0.0, k, th)
    q = 2 * k * np.sin(th / 2)
    a2.plot(np.degrees(th), np.abs(f) ** 2 * 10, color=c, ls="-")
    a2.plot(np.degrees(th), fb(V0, q) ** 2 * 10, color=c, ls="--")
a2.set(yscale="log", xlabel=r"$\theta_\mathrm{cm}$ [graus]", ylabel=r"$\mathrm{d}\sigma/\mathrm{d}\Omega$ [mb/sr]", xlim=(0, 180))
a2.set_xticks([0, 45, 90, 135, 180])
a2.text(0.97, 0.95, f"$E_\\mathrm{{cm}} = {Ed:g}$ MeV", transform=a2.transAxes, ha="right", va="top", fontsize=8)
virgula(a2, "x")
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout(rect=(0, 0.12, 1, 1))
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
salvar(fig, "5.3/fig02_validade_born_gaussiano_v2")
