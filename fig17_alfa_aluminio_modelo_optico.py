"""Figura 14 do TCC (Seção 4.4.3) – Modelo óptico completo (Coulomb + Woods-Saxon complexo) para alfa + 27Al
comparado aos dados de Bieler (1924) (Tabela do Apêndice D) e predição da função de excitação a ângulo fixo
(arranjo de Rutherford & Chadwick, 1925).

Potencial: U(r) = V_C(r; R_C) - V f(r;R,a) - i W f(r;R,a),  R = 1,40 A^{1/3}, a = 0,52 fm,
V = 185 MeV, W = 25 MeV (potencial de McFadden & Satchler, 1966), R_C = 1,3 A^{1/3} fm.
Coulomb: funções F_l, G_l (Steed); S_l nuclear por casamento em r > alcance nuclear.
Os dados de Bieler são razões de contagens integradas sobre os diafragmas anulares; o modelo é
integrado sobre a mesma faixa angular (laboratório) para comparação justa.
"""
import numpy as np
from estilo import *
from scatter_core import *

A, Z = 27, 13
mT = A * AMU; mu = M_ALPHA * mT / (M_ALPHA + mT)
V, W, r0, a, rc = 185.0, 25.0, 1.40, 0.52, 1.3

def calc(Elab, V=V, W=W, lmax=30):
    Ecm = Elab * mT / (M_ALPHA + mT)
    k = np.sqrt(2 * mu * Ecm) / HBARC
    eta = 2 * Z * E2 * mu / (HBARC**2 * k)
    R = r0 * A ** (1 / 3); Rc = rc * A ** (1 / 3)
    U = lambda r: coulomb_sphere(r, 2 * Z, Rc) - (V + 1j * W) * woods_saxon(r, R, a)
    rm = max(R + 14 * a, 2.5 * eta / k + 2)
    S, k = smatrix(U, mu, Ecm, lmax, rm, eta=eta, h=0.01)
    return S, k, eta

def ratio_lab(S, k, eta, th_lab):
    return ratio_to_rutherford(S, eta, k, lab_to_cm_angle(th_lab, M_ALPHA, mT))

bins = np.array([[17.6, 23.6], [23.7, 44.2], [43.9, 59.7], [60.3, 99.5]])
th_ef = np.array([20.2, 30.9, 50.7, 75.8])  # média ponderada por cossec^3(theta/2) (Bieler 1924, Eq. 1)
dat = np.array([1.00, 0.94, 0.78, 0.71]); err = np.array([0.03, 0.02, 0.02, 0.02])

def binned(S, k, eta):
    """Razão sigma/sigma_R média em cada diafragma, ponderada pelo número de partículas que a lei de Rutherford
    prevê na geometria de Bieler (folha anular, fonte e anteparo no eixo): peso cossec^3(theta_lab/2)
    (Bieler 1924, p. 438, Eq. 1; Apêndice D). É a mesma média usada na Figura 9."""
    out = []
    for b in np.radians(bins):
        tl = np.linspace(*b, 600); w = np.sin(tl / 2) ** -3
        out.append(np.trapezoid(ratio_lab(S, k, eta, tl) * w, tl) / np.trapezoid(w, tl))
    return np.array(out)

if __name__ == "__main__":
    # convergência em lmax e no raio de casamento
    S1, k, eta = calc(7.68, lmax=25); S2, _, _ = calc(7.68, lmax=40)
    th = np.radians(np.linspace(5, 175, 341))
    print("convergência lmax 25->40: max dif =", np.max(np.abs(ratio_lab(S1, k, eta, th) - ratio_lab(S2, k, eta, th))))
    print(f"eta = {eta:.3f}, k = {k:.3f} fm^-1;  |S_l| (l=0..8):", np.round(np.abs(S1[:9]), 3))

    # energias pelos alcances médios de Bieler (1924, Tab. II): E = E(7,0 cm) (R/7,0)^(2/3), E(7,0 cm) = 7,644 MeV (Apêndice D)
    Eal = lambda R: 7.644 * (R / 7.0) ** (2 / 3)
    E_RaC, E_Po = Eal(6.6), Eal(3.4)
    print(f"E(RaC, 6,6 cm) = {E_RaC:.3f} MeV;  E(Po, 3,4 cm) = {E_Po:.3f} MeV")
    fmtE = lambda E: f"{E:.2f}".replace(".", "{,}")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.95), gridspec_kw={"width_ratios": [1.35, 1]})
    thl = np.linspace(5, 172, 400)
    # faixa de energia (incerteza da energia efetiva no alvo)
    band = [ratio_lab(*calc(E), np.radians(thl)) for E in np.linspace(7.0, 7.68, 5)]
    a1.fill_between(thl, np.min(band, 0), np.max(band, 0), color=LARANJA, alpha=0.18, lw=0,
                    label=r"complexo, $E_\alpha = 7{,}0$–$7{,}68$ MeV")
    Sc, k, eta = calc(E_RaC)
    print(f"7,35 MeV: eta = {eta:.3f}, k = {k:.3f} fm^-1;  |S_l| (l=0..8):", np.round(np.abs(Sc[:9]), 3))
    a1.plot(thl, ratio_lab(Sc, k, eta, np.radians(thl)), color=LARANJA,
            label=r"complexo ($W=25$ MeV), $" + fmtE(E_RaC) + r"$ MeV")
    Sr, kr, er = calc(E_RaC, W=0.0)
    rr = ratio_lab(Sr, kr, er, np.radians(thl))
    print(f"real em 7,35 MeV: min = {rr.min():.3f} em {thl[np.argmin(rr)]:.0f} graus; max = {rr.max():.2f}")
    rcx = ratio_lab(Sc, k, eta, np.radians(thl)); imx = np.argmax(rcx)
    print(f"real em 7,35 MeV: maior pico abaixo de 165 graus = {rr[thl < 165].max():.2f}; o máximo de {rr.max():.2f} é a borda da malha ({thl[-1]:.0f} graus)")
    print(f"complexo em 7,35 MeV: máximo = {rcx[imx]:.3f} em {thl[imx]:.0f} graus; decrescente depois disso: {bool(np.all(np.diff(rcx[imx:]) < 0))}")
    a1.plot(thl, rr, color=AZUL, ls="--", label=r"real ($W=0$), $" + fmtE(E_RaC) + r"$ MeV")
    Sp, kp, ep = calc(E_Po)
    a1.plot(thl, ratio_lab(Sp, kp, ep, np.radians(thl)), color=VERDE, lw=1.2,
            label=r"complexo, $" + fmtE(E_Po) + r"$ MeV (polônio)")
    a1.axhline(1, color=CINZA, lw=0.8, ls=":")
    xerr = np.abs(bins.T - th_ef)
    a1.errorbar(th_ef, dat, yerr=err, xerr=xerr, fmt="o", color=TINTA, ms=5, capsize=2, lw=0.9,
                label=r"Bieler (1924), RaC (alcance 6,6 cm)", zorder=8)
    datp, errp = np.array([1.00, 0.93]), np.array([0.03, 0.03])
    a1.errorbar(th_ef[2:], datp, yerr=errp, xerr=xerr[:, 2:], fmt="D", color=TINTA, mfc="white", ms=4.5, capsize=2,
                lw=0.9, label=r"Bieler (1924), Po (alcance 3,4 cm)", zorder=8)
    mb = binned(Sc, k, eta)
    a1.plot(th_ef, mb, "s", mfc="white", mec=LARANJA, ms=6, mew=1.3, label="modelo integrado nos diafragmas", zorder=6)
    mbp = binned(Sp, kp, ep)[2:]
    a1.plot(th_ef[2:], mbp, "s", mfc="white", mec=VERDE, ms=6, mew=1.3, zorder=6)
    chi2 = np.sum(((mb - dat) / err) ** 2)
    print("Bieler RaC, modelo integrado (7,35 MeV):", np.round(mb, 3), " chi2 =", round(chi2, 1))
    print("Bieler Po, modelo integrado (4,72 MeV), III e IV:", np.round(mbp, 3))
    for E in [7.0, 7.68]:
        m = binned(*calc(E)); print(f"   E={E}: binned", np.round(m, 3), " chi2 =", round(np.sum(((m - dat) / err) ** 2), 1))
    lr = np.arange(len(Sc)); srx = np.pi / k**2 * np.sum((2 * lr + 1) * (1 - np.abs(Sc) ** 2)) * 10
    print(f"sigma_reacao(7,35) = {srx:.0f} mb")
    a1.set(xlabel=r"$\theta_\mathrm{lab}$ [graus]", ylabel=r"$\sigma/\sigma_\mathrm{R}$", xlim=(0, 175), ylim=(0, 1.6))
    a1.set_xticks([0, 30, 60, 90, 120, 150])
    h1, l1 = a1.get_legend_handles_labels()

    # função de excitação a ângulo fixo (Rutherford & Chadwick 1925)
    Es = np.linspace(3.0, 9.0, 61)
    for thf, c, ls in [(135, TINTA, "-"), (90, TINTA, "-.")]:
        r = [ratio_lab(*calc(E), np.radians([thf]))[0] for E in Es]
        a2.plot(Es, r, color=c, ls=ls, label=fr"$\theta_\mathrm{{lab}} = {thf}^\circ$")
    Eb = 2 * Z * E2 / (1.2 * (4 ** (1 / 3) + A ** (1 / 3))) * (M_ALPHA + mT) / mT
    a2.axvline(Eb, color=CINZA, lw=0.8, ls=":")
    a2.text(Eb - 0.1, 1.52, "barreira no raio\nde toque:\n" + (f"${Eb:.1f}$ MeV").replace(".", "{,}"),
            fontsize=7, color=CINZA, ha="right", va="top")
    print(f"barreira de toque: E_lab = {Eb:.2f} MeV")
    # barreira l = 0 do próprio potencial óptico (Coulomb + parte real nuclear)
    from scipy.optimize import minimize_scalar
    R0, Rc0 = r0 * A ** (1 / 3), rc * A ** (1 / 3)
    Vef = lambda r: coulomb_sphere(r, 2 * Z, Rc0) - V * woods_saxon(r, R0, a)
    opt = minimize_scalar(lambda r: -Vef(r), bounds=(5, 12), method="bounded")
    EbO = Vef(opt.x) * (M_ALPHA + mT) / mT
    print(f"barreira l=0 do potencial óptico: r_B = {opt.x:.2f} fm, V_B(cm) = {Vef(opt.x):.2f} MeV, E_lab = {EbO:.2f} MeV")
    a2.axvline(EbO, color=CINZA, lw=0.8, ls="--")
    a2.text(EbO - 0.1, 1.52, "barreira $l=0$ do\npotencial óptico:\n" + (f"${EbO:.1f}$ MeV").replace(".", "{,}"),
            fontsize=7, color=CINZA, ha="right", va="top")
    a2.axhline(1, color=CINZA, lw=0.8, ls=":")
    a2.set(xlabel=r"$E_\alpha$ (lab) [MeV]", ylabel=r"$\sigma/\sigma_\mathrm{R}$", ylim=(0, 1.6), xlim=(3, 9))
    a2.legend(fontsize=7.5, loc="lower left")
    for ax in (a1, a2): virgula(ax)
    rotulo(a1, "a"); rotulo(a2, "b")
    fig.tight_layout(rect=(0, 0.175, 1, 1))
    fig.legend(h1, l1, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=7, title_fontsize=7, frameon=False,
               title=r"Painel (a) $-$ potencial de Woods--Saxon de volume: $V = 185$ MeV, $W = 25$ MeV, "
                     r"$r_0 = 1{,}40$ fm, $a = 0{,}52$ fm; $r_C = 1{,}3$ fm")
    salvar(fig, "6.4/fig08_alfa_al_modelo_optico_bieler_v3")
    i = np.argmin(np.abs(Es - 5.0)); print("135 graus: razão em 5 MeV e 7,68 MeV:",
          round(ratio_lab(*calc(5.0), np.radians([135]))[0], 3), round(ratio_lab(*calc(7.68), np.radians([135]))[0], 3))
