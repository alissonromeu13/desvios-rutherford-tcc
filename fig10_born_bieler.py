"""Figura 10 do TCC (Seção 4.1.2) – Interferência Coulomb-nuclear na 1ª aproximação de Born (poço quadrado real).

sigma/sigma_R = [1 - g x j1(x)]^2,  x = qR = 2kR sen(theta_cm/2),  g = V0/V_C(R) = V0 R/(k_e z Z e^2).

(a) Modelo puro: kR = 3,5, K = 0,35 (em função de theta_cm).
(b) Dados de Bieler (1924), Tabs. II (Al, RaC, alcance médio 6,6 cm) e III (Mg, RaC, 6,2 cm), p. 445;
    faixas angulares da Tab. I (p. 445). Cada ponto é uma média sobre uma faixa de ângulo de laboratório:
    o modelo é comparado com a mesma média, ponderada pelo número de partículas de Rutherford na geometria
    de Bieler (folha anular, fonte e anteparo no eixo): peso cossec^3(theta/2) (Bieler 1924, p. 438, Eq. 1; Apêndice F). Conversão lab -> CM exata para colisão elástica (só depende de m/M).
    kR_ef e K são parâmetros de AJUSTE (mínimos quadrados ponderados), não grandezas derivadas.
    Al: ajuste de dois parâmetros. Mg (3 pontos): kR_ef = (kR_ef/kR_fís do Al) × kR_fís do Mg, ajusta-se só K.
"""
import numpy as np
from scipy.special import spherical_jn
from scipy.optimize import least_squares, brentq
from scipy.integrate import quad
from estilo import *

M_ALFA = 4.002602
ALVOS = {  # massa (u), Z, dados: (theta_min, theta_max, razão, erro) — Bieler (1924), Tabs. I–III, p. 445
    "Al": dict(M=26.981538, Z=13, alcance=6.6, dados=[(17.6, 23.6, 1.00, 0.03), (23.7, 44.2, 0.94, 0.02),
                                                     (43.9, 59.7, 0.78, 0.02), (60.3, 99.5, 0.71, 0.02)]),
    "Mg": dict(M=23.985042, Z=12, alcance=6.2, dados=[(17.6, 23.6, 1.02, 0.03),
                                                     (43.9, 59.7, 0.90, 0.03), (60.3, 99.5, 0.69, 0.02)]),
}

def razao(theta_cm, kR, K):
    x = 2 * kR * np.sin(np.asarray(theta_cm) / 2)
    return (1 - K * x * spherical_jn(1, x)) ** 2

def lab2cm(th_lab, tau):
    """theta_cm a partir de theta_lab (rad), tau = m/M < 1 (ramo único)."""
    return th_lab + np.arcsin(tau * np.sin(th_lab))

def jacobiano(th_cm, tau):
    """dOmega_cm/dOmega_lab."""
    return (1 + tau**2 + 2 * tau * np.cos(th_cm)) ** 1.5 / abs(1 + tau * np.cos(th_cm))

def peso_rutherford_lab(th_lab, tau):
    """Peso de Rutherford na geometria de Bieler (p. 438, Eq. 1, fórmula de Chadwick): folha anular com
    fonte e anteparo no eixo, a distâncias iguais. N(phi1, phi2) ∝ G(phi1) − G(phi2), com
    G(phi) = cot(phi/2) cosec(phi/2) − ln tan(phi/4); o peso por unidade de ângulo é −dG/dphi = cosec^3(phi/2).
    Bieler despreza o recuo nessa fórmula (correção < 2 % até 100°, p. 439), então o peso é o de laboratório
    sem recuo; tau entra apenas na conversão do ângulo do modelo."""
    p = th_lab
    return np.sin(p / 2) ** -3

def media_faixa(f, a, b, tau):
    a, b = np.radians(a), np.radians(b)
    num = quad(lambda t: f(lab2cm(t, tau)) * peso_rutherford_lab(t, tau), a, b, limit=200)[0]
    den = quad(lambda t: peso_rutherford_lab(t, tau), a, b, limit=200)[0]
    return num / den

def angulo_efetivo(a, b, tau):
    a, b = np.radians(a), np.radians(b)
    num = quad(lambda t: t * peso_rutherford_lab(t, tau), a, b)[0]
    return np.degrees(num / quad(lambda t: peso_rutherford_lab(t, tau), a, b)[0])

def k_fisico(M, alcance):
    """Energia pela lei de Geiger (v^3 ∝ alcance), com a mesma âncora do Apêndice F e da Figura 17:
    v0 = 1,92e9 cm/s para 7,0 cm de ar, isto é, E(7,0 cm) = 7,644 MeV (Evans 1955, p. 512, Tab. 1.1, arredondado)."""
    E = 7.644 * (alcance / 7.0) ** (2 / 3)
    mu = M_ALFA * M / (M_ALFA + M) * 931.494
    Ecm = E * M / (M_ALFA + M)
    return E, Ecm, np.sqrt(2 * mu * Ecm) / 197.327

res = {}
for nome, d in ALVOS.items():
    tau = M_ALFA / d["M"]
    dados = np.array(d["dados"])
    def resid(p):
        return [(media_faixa(lambda t: razao(t, *p), a, b, tau) - y) / s for a, b, y, s in dados]
    E, Ecm, k = k_fisico(d["M"], d["alcance"])
    R = 1.2 * (d["M"] ** (1 / 3))
    if nome == "Al":   # ajuste de dois parâmetros (kR_ef e K)
        fit = least_squares(resid, x0=[1.7, 0.17], bounds=([0.3, 0.0], [6.0, 2.0]))
        J = fit.jac; cov = np.linalg.inv(J.T @ J); ndf = len(dados) - 2
        fator = fit.x[0] / (k * R)   # kR_ef / kR_fisico, reaproveitado no Mg
    else:              # Mg: três pontos não determinam dois parâmetros; kR_ef = mesmo fator × kR_fisico
        kRf = fator * k * R
        r1 = lambda p: resid([kRf, p[0]])
        f1 = least_squares(r1, x0=[0.2], bounds=([0.0], [5.0]))
        fit = type("F", (), {})(); fit.x = np.array([kRf, f1.x[0]]); fit.fun = f1.fun
        c = np.linalg.inv(f1.jac.T @ f1.jac)[0, 0]; cov = np.array([[0, 0], [0, c]]); ndf = len(dados) - 1
    chi2 = np.sum(fit.fun**2)
    res[nome] = dict(tau=tau, p=fit.x, err=np.sqrt(np.diag(cov)), chi2=chi2, ndf=ndf, dados=dados)
    print(f"{nome}: kR_ef = {fit.x[0]:.3f} ± {np.sqrt(cov[0,0]):.3f}, K = {fit.x[1]:.3f} ± {np.sqrt(cov[1,1]):.3f}, "
          f"chi2/ndf = {chi2:.2f}/{ndf}, K kR^2 = {fit.x[1]*fit.x[0]**2:.3f}")
    if nome == "Al":   # robustez do ajuste de dois parâmetros: outros pontos de partida
        alt = [least_squares(resid, x0=x0, bounds=([0.3, 0.0], [6.0, 2.0])).x for x0 in ([1.0, 0.4], [3.0, 0.1], [2.5, 0.3])]
        print("   outros pontos de partida ->", [tuple(float(v) for v in np.round(a, 3)) for a in alt])
    print(f"   E_lab ≈ {E:.2f} MeV, E_cm ≈ {Ecm:.2f} MeV, k = {k:.3f} fm^-1, kR(R=1,2A^1/3={R:.2f} fm) = {k*R:.2f}")
    for a, b, y, s in dados:
        print(f"   faixa {a}–{b}°: theta_ef(lab) = {angulo_efetivo(a, b, tau):.1f}°, dado {y} ± {s}, "
              f"modelo = {media_faixa(lambda t: razao(t, *fit.x), a, b, tau):.3f}")
    xmax = 2.7437
    thcm_min = 2 * np.degrees(np.arcsin(min(1, xmax / (2 * fit.x[0]))))
    print(f"   mínimo do modelo em theta_cm = {thcm_min:.0f}°")
    if nome == "Al":   # números citados no texto da Seção 4.1.2
        r1 = lambda p: resid([k * R, p[0]])               # kR fixo no valor físico; ajusta-se só g
        f1 = least_squares(r1, x0=[0.2], bounds=([0.0], [5.0]))
        g1 = f1.x[0]; x1 = brentq(lambda x: x * spherical_jn(1, x), 3.0, 5.0)   # razão volta a 1 no zero de j1
        print(f"   kR físico fixo = {k*R:.2f}: g = {g1:.3f}, chi2 = {np.sum(f1.fun**2):.0f} para {len(dados)-1} graus de liberdade;"
              f" razão volta a 1 em theta_cm = {2*np.degrees(np.arcsin(x1/(2*k*R))):.0f}°")
        print(f"   80° lab -> {np.degrees(lab2cm(np.radians(80.0), tau)):.1f}° cm;"
              f"  V_C(4 fm) = {2*d['Z']*1.43996/4:.1f} MeV;  k_e zZ e^2/R = {2*d['Z']*1.43996/R:.1f} MeV")
        Rf = fit.x[0] / k
        print(f"   V0 = g k_e zZ e^2/R: {fit.x[1]*2*d['Z']*1.43996/R:.1f} MeV (R físico),"
              f" {fit.x[1]*2*d['Z']*1.43996/Rf:.1f} MeV (R_ef = kR_ef/k = {Rf:.2f} fm)")

# ---------------- figura ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3))
th = np.linspace(0.5, 120, 1200)
a1.axhline(1, color=CINZA, lw=0.9, ls="--", label="Rutherford")
a1.plot(th, razao(np.radians(th), 3.5, 0.35), color=AZUL, label=r"Born, $kR = 3{,}5$, $g = 0{,}35$")
thz = 2 * np.degrees(np.arcsin(4.4934 / 7.0)); thm = 2 * np.degrees(np.arcsin(2.7437 / 7.0))
a1.axvline(thz, color=CINZA, lw=0.6, ls=(0, (1, 2)))
a1.text(thz + 1.5, 0.27, r"$qR = 4{,}49$" + "\n" + r"(zero de $j_1$)", fontsize=7.5, color=CINZA, va="bottom")
a1.plot([thm], [razao(np.radians(thm), 3.5, 0.35)], "o", ms=3.5, color=AZUL)
a1.text(thm, razao(np.radians(thm), 3.5, 0.35) - 0.07, r"$qR = 2{,}74$", fontsize=7.5, color=AZUL, ha="center", va="top")
a1.set(xlim=(0, 120), ylim=(0.2, 1.6), xlabel=r"$\theta_{\text{cm}}$ (graus)", ylabel=r"$\sigma/\sigma_{\text{R}}$")
a1.legend(loc="upper left", fontsize=7.5)

thl = np.linspace(1, 140, 800)
a2.axhline(1, color=CINZA, lw=0.9, ls="--")
estilo_alvo = {"Al": (AZUL, "o", "-"), "Mg": (LARANJA, "s", "--")}
for nome, r in res.items():
    cor, mk, ls = estilo_alvo[nome]
    tau = r["tau"]; kR, K = r["p"]
    a2.plot(thl, razao(lab2cm(np.radians(thl), tau), kR, K), color=cor, ls=ls, lw=1.3,
            label=(r"$^{27}$Al" if nome == "Al" else r"$^{24}$Mg") + rf": $kR_{{\text{{ef}}}} = {kR:.2f}$, $g = {K:.2f}$".replace(".", "{,}"))
    for a, b, y, s in r["dados"]:
        te = angulo_efetivo(a, b, tau)
        a2.errorbar(te, y, yerr=s, xerr=[[te - a], [b - te]], fmt=mk, ms=4, color=cor, mfc=cor,
                    elinewidth=0.8, capsize=1.5)
        a2.plot(te, media_faixa(lambda t: razao(t, kR, K), a, b, tau), mk, ms=6, mfc="none", mec=cor, mew=0.9)
a2.plot([], [], "o", color=TINTA, ms=4, label="Bieler (1924), média na faixa")
a2.plot([], [], "o", mfc="none", mec=TINTA, ms=6, label="modelo, média na faixa")
a2.set(xlim=(0, 140), ylim=(0.5, 1.12), xlabel=r"$\theta_{\text{lab}}$ (graus)", ylabel=r"$\sigma/\sigma_{\text{R}}$")
a2.legend(loc="lower left", fontsize=7)
for ax in (a1, a2):
    virgula(ax)
rotulo(a1, "a"); rotulo(a2, "b")
fig.tight_layout()
salvar(fig, "Potencial_combinado_e_superposicao_de_amplitudes/fig10_born_bieler")
