"""Núcleo numérico para as figuras do TCC (espalhamento não relativístico).

Unidades: MeV, fm, mb. Massas em MeV/c^2. Constantes do CODATA 2018 (Tiesinga et al. 2021): hbar*c = 197,327 MeV fm,
alfa = 1/137,036, m_alfa = 3727,379 MeV, m_n = 939,565 MeV; unidade de massa atômica 931,494 MeV.

Conteúdo
  coulomb_FG(eta, rho, lmax)   funções de Coulomb F_l, G_l e derivadas (método de Steed,
                               Barnett et al. 1974) -- testadas contra Bessel (eta=0) e Wronskiano.
  numerov_u(...)               integra u_l'' = [l(l+1)/r^2 + 2mu U(r)/hbar^2 - k^2] u_l
                               para todos os l ao mesmo tempo (vetorizado).
  smatrix(...)                 S_l nuclear por casamento com H^(+-) = G +- iF (ou Riccati-Hankel se eta=0).
  amplitudes(...)              f_C(theta), f_N(theta) e razão sigma/sigma_R.
"""
import numpy as np
from scipy.special import loggamma, spherical_jn, spherical_yn, eval_legendre

HBARC = 197.3269804
ALPHA = 1 / 137.035999
E2 = HBARC * ALPHA          # e^2/(4 pi eps0) = 1.44 MeV fm
AMU = 931.494
M_ALPHA = 3727.379
M_N = 939.565


# ---------------------------------------------------------------- Coulomb functions
def coulomb_FG(eta, rho, lmax, acc=1e-14):
    """F_l, F'_l, G_l, G'_l (derivadas em rho) para l = 0..lmax, rho > 0 real.

    CF1 (F'/F em lmax) + recorrência descendente; CF2 (p+iq) em l=0;
    normalização pelo Wronskiano F'G - FG' = 1; G sobe por recorrência.
    Válido para rho acima de ~ eta + sqrt(eta^2) (CF2 converge rápido aí).
    """
    L = lmax
    S = lambda l: l / rho + eta / l
    R = lambda l: np.sqrt(1 + (eta / l) ** 2)
    # CF1: f_L = S_{L+1} - R_{L+1}^2/(T_{L+1} - R_{L+2}^2/(T_{L+2} - ...)),  T_j = S_j + S_{j+1}
    # avaliado pelo método de Lentz modificado
    tiny = 1e-300
    f = S(L + 1)
    C, D = f, 0.0
    j = L + 1
    while True:
        a = -R(j) ** 2
        b = S(j) + S(j + 1)
        D = b + a * D
        D = tiny if D == 0 else 1 / D
        C = b + a / C
        C = tiny if C == 0 else C
        delta = C * D
        if j == L + 1:
            # primeiro passo: f = S(L+1) + a/(b + ...)
            C, D = b, 0.0
            tail, j = b, j + 1
            continue
        tail *= delta
        j += 1
        if abs(delta - 1) < acc or j > L + 200000:
            break
    f = S(L + 1) - R(L + 1) ** 2 / tail
    # downward recurrence (unnormalised)
    F = np.zeros(L + 1); Fp = np.zeros(L + 1)
    F[L] = 1.0
    Fp[L] = f * F[L]
    for l in range(L, 0, -1):
        F[l - 1] = (S(l) * F[l] + Fp[l]) / R(l)
        Fp[l - 1] = S(l) * F[l - 1] - R(l) * F[l]
        if abs(F[l - 1]) > 1e250:
            F[l - 1:] *= 1e-250; Fp[l - 1:] *= 1e-250
    f0 = Fp[0] / F[0]
    # CF2 (Steed) at l = 0
    xi = 1 / rho
    pq = complex(0.0, 1 - eta * xi)
    ar, ai = -eta * eta, eta          # a*b = -(eta^2 + l(l+1)) + i eta at l=0
    br, bi = 2 * (rho - eta), 2.0
    wi = 2 * eta
    dr = br / (br * br + bi * bi); di = -bi / (br * br + bi * bi)
    dp = -xi * (ar * di + ai * dr); dq = xi * (ar * dr - ai * di)
    p, q = 0.0, 1 - eta * xi
    pk = 0.0
    for _ in range(200000):
        p += dp; q += dq
        pk += 2
        ar += pk; ai += wi; bi += 2
        d = ar * dr - ai * di + br
        di = ai * dr + ar * di + bi
        c = 1 / (d * d + di * di)
        dr = c * d; di = -c * di
        a = br * dr - bi * di - 1
        b = bi * dr + br * di
        c = dp * a - dq * b
        dq = dp * b + dq * a
        dp = c
        if abs(dp) + abs(dq) < (abs(p) + abs(q)) * acc:
            break
    else:
        raise RuntimeError("CF2 no convergence (rho too small?)")
    # normalisation: F0^2 = q / ((f-p)^2 + q^2)
    F0 = np.sign(F[0]) * np.sqrt(q / ((f0 - p) ** 2 + q * q))
    scale = F0 / F[0]
    F *= scale; Fp *= scale
    G = np.zeros(L + 1); Gp = np.zeros(L + 1)
    G[0] = F[0] * (f0 - p) / q
    Gp[0] = p * G[0] - q * F[0]
    for l in range(0, L):
        G[l + 1] = (S(l + 1) * G[l] - Gp[l]) / R(l + 1)
        Gp[l + 1] = R(l + 1) * G[l] - S(l + 1) * G[l + 1]
    return F, Fp, G, Gp


def sigma_l(eta, lmax):
    """Fases coulombianas sigma_l = arg Gamma(l+1+i eta)."""
    l = np.arange(lmax + 1)
    return np.imag(loggamma(l + 1 + 1j * eta))


# ---------------------------------------------------------------- Numerov
def numerov_u(Ufun, mu, E, lmax, rmax, h=0.01, rmin=None):
    """u_l(r) em toda a malha para l=0..lmax (vetorizado em l).

    Ufun(r): potencial total (complexo permitido) em MeV, incluindo Coulomb se houver.
    Retorna r (N,), u (N, lmax+1) complexos, normalização arbitrária por l.
    """
    k2 = 2 * mu * E / HBARC ** 2
    r = np.arange(1, int(round(rmax / h)) + 1) * h   # r_1 = h (u(0)=0)
    l = np.arange(lmax + 1)
    g = k2 - 2 * mu * Ufun(r)[:, None] / HBARC ** 2 - l[None, :] * (l[None, :] + 1) / r[:, None] ** 2
    w = 1 + h * h * g / 12
    u = np.zeros((len(r), lmax + 1), dtype=complex)
    # start: u ~ r^{l+1} at the first two points (Numerov needs u_0 = 0 at r=0)
    u[0] = h ** (l + 1) * 1.0
    u[1] = (2 * h) ** (l + 1) * 1.0
    # proper start with u(0)=0: use series for first step from r0=0: u1 given, u2 from Numerov with u0=0 (only exact for l>=1 when w0 finite);
    for i in range(1, len(r) - 1):
        u[i + 1] = ((12 - 10 * w[i]) * u[i] - w[i - 1] * u[i - 1]) / w[i + 1]
        big = np.abs(u[i + 1]) > 1e200
        if big.any():
            u[: i + 2, big] *= 1e-200
    return r, u


def smatrix(Ufun, mu, E, lmax, rmatch, eta=0.0, h=0.01):
    """S_l nuclear: u_l = A [H^-_l(kr) - S_l H^+_l(kr)], casado em dois pontos r1<r2 ~ rmatch.

    Para eta = 0 usa Riccati-Hankel; senão funções de Coulomb (Steed).
    Ufun deve conter a parte coulombiana quando eta != 0 e coincidir com o Coulomb
    pontual para r >= rmatch.
    """
    k = np.sqrt(2 * mu * E) / HBARC
    r, u = numerov_u(Ufun, mu, E, lmax, rmatch + 3 * h, h=h)
    i2 = len(r) - 1
    i1 = i2 - int(round(1.0 / h))   # 1 fm de separação
    r1, r2 = r[i1], r[i2]
    def H(rr):
        rho = k * rr
        if eta == 0.0:
            ll = np.arange(lmax + 1)
            F = rho * spherical_jn(ll, rho); G = -rho * spherical_yn(ll, rho)
        else:
            F, _, G, _ = coulomb_FG(eta, rho, lmax)
        return G - 1j * F, G + 1j * F    # H^-, H^+  (H^+ ~ e^{+i theta})
    Hm1, Hp1 = H(r1); Hm2, Hp2 = H(r2)
    ratio = u[i1] / u[i2]
    S = (Hm1 - ratio * Hm2) / (Hp1 - ratio * Hp2)
    return S, k


def coulomb_amplitude(eta, k, theta):
    s2 = np.sin(theta / 2) ** 2
    sig0 = sigma_l(eta, 0)[0]
    return -eta / (2 * k * s2) * np.exp(-1j * eta * np.log(s2) + 2j * sig0)


def nuclear_amplitude(S, eta, k, theta):
    l = np.arange(len(S))
    sl = sigma_l(eta, len(S) - 1) if eta != 0 else np.zeros(len(S))
    x = np.cos(theta)
    P = np.array([eval_legendre(li, x) for li in l])   # (L, Ntheta)
    coef = (2 * l + 1) * np.exp(2j * sl) * (S - 1) / (2j * k)
    return coef @ P


def ratio_to_rutherford(S, eta, k, theta):
    fC = coulomb_amplitude(eta, k, theta)
    fN = nuclear_amplitude(S, eta, k, theta)
    return np.abs(fC + fN) ** 2 / np.abs(fC) ** 2


# ---------------------------------------------------------------- potentials
def woods_saxon(r, R, a):
    return 1 / (1 + np.exp((r - R) / a))


def coulomb_sphere(r, zZ, Rc):
    r = np.asarray(r, float)
    return np.where(r < Rc, zZ * E2 / (2 * Rc) * (3 - (r / Rc) ** 2), zZ * E2 / np.maximum(r, 1e-12))


def lab_to_cm_angle(theta_lab, m1, m2):
    """theta_lab em rad -> theta_cm (elástico, m1 projétil)."""
    x = m1 / m2
    return theta_lab + np.arcsin(x * np.sin(theta_lab))


# ---------------------------------------------------------------- self tests
if __name__ == "__main__":
    # 1) eta = 0 must reproduce Riccati-Bessel
    rho = 7.3
    F, Fp, G, Gp = coulomb_FG(0.0, rho, 10)
    ll = np.arange(11)
    Fb = rho * spherical_jn(ll, rho); Gb = -rho * spherical_yn(ll, rho)
    print("eta=0  max|F-Fbessel| =", np.max(np.abs(F - Fb)), " max|G-Gbessel| =", np.max(np.abs(G - Gb)))
    # 2) Wronskian and known value check (eta=3, rho=15)
    for eta, rho in [(3.0, 15.0), (2.96, 12.7), (18.0, 60.0), (0.5, 3.0)]:
        F, Fp, G, Gp = coulomb_FG(eta, rho, 30)
        W = Fp * G - F * Gp
        # independent check: integrate pure Coulomb eq. for l=0 between rho and rho+1 with Numerov from (F,F') at rho
        print(f"eta={eta} rho={rho}: max|W-1| = {np.max(np.abs(W-1)):.2e}")
    # 3) numerical vs analytic F_l: derivative consistency F'(rho) by finite difference
    eta, rho, d = 3.0, 15.0, 1e-5
    Fa, Fpa, Ga, Gpa = coulomb_FG(eta, rho, 5)
    Fb_, _, Gb_, _ = coulomb_FG(eta, rho + d, 5); Fc_, _, Gc_, _ = coulomb_FG(eta, rho - d, 5)
    print("finite-diff F' error:", np.max(np.abs((Fb_ - Fc_) / (2 * d) - Fpa)), " G' error:", np.max(np.abs((Gb_ - Gc_) / (2 * d) - Gpa)))
    # 4) pure point Coulomb (R_c tiny) -> S_l = 1
    mu = M_ALPHA * 27 * AMU / (M_ALPHA + 27 * AMU); Ecm = 6.69
    k = np.sqrt(2 * mu * Ecm) / HBARC; eta = 2 * 13 * E2 * mu / (HBARC ** 2 * k)
    S, k = smatrix(lambda r: coulomb_sphere(r, 26, 0.05), mu, Ecm, 20, 30.0, eta=eta, h=0.005)
    print("pure Coulomb: max|S_l - 1| =", np.max(np.abs(S - 1)), " eta=", eta)
    # 5) neutral square-well s-wave vs analytic
    V0, Rw = 40.0, 3.0; mu = M_N; E = 5.0
    S, k = smatrix(lambda r: np.where(r < Rw, -V0, 0.0) + 0j, mu, E, 0, 10.0, eta=0.0, h=0.001)
    K = np.sqrt(2 * mu * (E + V0)) / HBARC
    d0 = -k * Rw + np.arctan(k / K * np.tan(K * Rw))
    print("square well delta0 numeric =", np.angle(S[0]) / 2, " analytic =", (d0 + np.pi / 2) % np.pi - np.pi / 2)
