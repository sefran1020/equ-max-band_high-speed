"""
jitter.py — Apertura HORIZONTAL del ojo: curva de bañera y jitter (RJ + DJ).

Fase α del plan `planificacion/plan_jitter_y_coopt.md`. Reutiliza el pipeline
validado (`cadena_enlace`) sin modificarlo. La idea, semi-analítica y coherente
con el BER vertical existente:

  1) Curva de bañera CLEAN: BER en función de la fase de muestreo dentro de la UI
     (se muestrea la MISMA forma de onda en cada fase y se reusa el BER por niveles).
  2) Jitter del reloj de muestreo: la BER observada a una fase nominal es la
     esperanza sobre el jitter -> convolución (circular, el ojo es periódico en UI)
     de la bañera con la PDF de jitter (dual-Dirac: RJ gaussiano + DJ pico-a-pico).
  3) Apertura horizontal (eye width) = tramo de UI con BER <= umbral.

No incluye DFE: la bañera se evalúa sobre la forma de onda ecualizada linealmente
(entrada del decisor), que es donde se mide la apertura horizontal.
"""

import numpy as np

from cadena_enlace import (
    SPS, Q, estadisticas_niveles, _decidir_pam4,
)


# --------------------------------------------------------------------------- #
# Alineación: offset absoluto del muestreo óptimo del símbolo 0
# --------------------------------------------------------------------------- #
def alinear_base(y, sym, sps_val=SPS, n_score=400):
    """Devuelve el offset de muestra (absoluto) del muestreo óptimo del símbolo 0,
    es decir base = fase_optima + retardo_en_simbolos * sps."""
    n = len(y)
    sym = np.asarray(sym, dtype=float)
    x_up = np.repeat(sym, sps_val)
    x_up = x_up[:n] if len(x_up) >= n else np.pad(x_up, (0, n - len(x_up)))
    cc = np.fft.irfft(np.fft.rfft(y - y.mean()) * np.conj(np.fft.rfft(x_up - x_up.mean())), n=n)
    off0 = int(np.argmax(cc[: n // 2]))
    lag0 = off0 // sps_val
    mejor = None
    for ph in range(sps_val):
        ds = y[ph::sps_val]
        for lag in (lag0 - 1, lag0, lag0 + 1):
            if lag < 0:
                continue
            navail = min(n_score, len(ds) - lag)   # ventana disponible tras el retardo
            if navail < 50:
                continue
            a = ds[lag:lag + navail]
            a = a - a.mean()
            b = sym[:navail] - sym[:navail].mean()
            na = np.linalg.norm(a)
            nb = np.linalg.norm(b)
            if na <= 0 or nb <= 0:
                continue
            score = abs(float(np.dot(a, b)) / (na * nb))
            if mejor is None or score > mejor[0]:
                mejor = (score, ph, lag)
    if mejor is None:
        return off0                                # fallback: offset de resolución completa
    _, ph, lag = mejor
    return ph + lag * sps_val


# --------------------------------------------------------------------------- #
# BER vertical a partir de muestras (sin DFE)
# --------------------------------------------------------------------------- #
def ber_pam4_muestras(m, tx):
    mu, s, cnt = estadisticas_niveles(m, tx)
    if np.any(cnt == 0):
        return 1.0
    o = np.argsort(mu)
    mu, s, cnt = mu[o], s[o], cnt[o]
    thr = (mu[:-1] + mu[1:]) / 2.0
    eps = 1e-300
    Pe = np.zeros(4)
    Pe[0] = Q((thr[0] - mu[0]) / (s[0] + eps))
    Pe[1] = Q((mu[1] - thr[0]) / (s[1] + eps)) + Q((thr[1] - mu[1]) / (s[1] + eps))
    Pe[2] = Q((mu[2] - thr[1]) / (s[2] + eps)) + Q((thr[2] - mu[2]) / (s[2] + eps))
    Pe[3] = Q((mu[3] - thr[2]) / (s[3] + eps))
    w = cnt / cnt.sum()
    return float(np.sum(w * Pe)) / 2.0


# --------------------------------------------------------------------------- #
# Curva de bañera CLEAN (BER vs fase) sobre una UI centrada en el óptimo
# --------------------------------------------------------------------------- #
def banera_clean(y, sym, sps_val=SPS, base=None, n_descartar=100):
    if base is None:
        base = alinear_base(y, sym, sps_val)
    n = len(y)
    sym = np.asarray(sym, dtype=float)
    Nsym = len(sym)
    deltas = np.arange(-sps_val // 2, sps_val // 2)   # muestras, centrado en 0
    bers = np.empty(len(deltas))
    for i, d in enumerate(deltas):
        off = base + d
        idx = off + sps_val * np.arange(Nsym)
        idx = idx[(idx >= 0) & (idx < n)]
        k = len(idx)
        if k <= n_descartar + 20:
            bers[i] = 1.0
            continue
        m = y[idx][n_descartar:]
        tx = sym[:k][n_descartar:]
        bers[i] = ber_pam4_muestras(m, tx)
    fase_ui = deltas / sps_val
    return fase_ui, np.clip(bers, 1e-12, 1.0)


def ber_pam4_dfe(m, tx, n_dfe):
    """BER vertical con DFE data-aided (mínimos cuadrados) aplicado a las muestras."""
    nt = n_dfe
    N = len(m)
    if N <= nt + 20:
        return 1.0
    Xreg = np.column_stack([tx[nt - i: N - i] for i in range(nt + 1)])
    h, *_ = np.linalg.lstsq(Xreg, m[nt:N], rcond=None)
    h0 = h[0] if abs(h[0]) > 1e-12 else 1.0
    post = h[1:]
    hist = [0.0] * nt
    m_eq = np.empty_like(m)
    for kk in range(N):
        fb = sum(post[i] * hist[i] for i in range(nt))
        yk = m[kk] - fb
        m_eq[kk] = yk
        hist = [_decidir_pam4(yk, h0)] + hist[:-1]
    return ber_pam4_muestras(m_eq, tx)


def banera_dfe(y, sym, sps_val=SPS, n_dfe=2, base=None, n_descartar=100):
    """Bañera DFE-AWARE: en cada fase se aplica el DFE data-aided antes del BER.
    Refleja el ojo del enlace COMPLETO (post-DFE), a diferencia de `banera_clean`
    (entrada del decisor lineal, pre-DFE)."""
    if base is None:
        base = alinear_base(y, sym, sps_val)
    n = len(y)
    sym = np.asarray(sym, dtype=float)
    Nsym = len(sym)
    deltas = np.arange(-sps_val // 2, sps_val // 2)
    bers = np.empty(len(deltas))
    for i, d in enumerate(deltas):
        off = base + d
        idx = off + sps_val * np.arange(Nsym)
        idx = idx[(idx >= 0) & (idx < n)]
        kk = len(idx)
        if kk <= n_descartar + 40:
            bers[i] = 1.0
            continue
        m = y[idx][n_descartar:]
        tx = sym[:kk][n_descartar:]
        bers[i] = ber_pam4_dfe(m, tx, n_dfe)
    return deltas / sps_val, np.clip(bers, 1e-12, 1.0)


def refinar(fase_ui, ber, n_fine=512):
    """Interpola log10(BER) de forma periódica (1 UI) a una rejilla fina."""
    lb = np.log10(np.clip(ber, 1e-12, 1.0))
    order = np.argsort(fase_ui)
    ph, lb = fase_ui[order], lb[order]
    ph_ext = np.concatenate([ph - 1.0, ph, ph + 1.0])
    lb_ext = np.concatenate([lb, lb, lb])
    xf = np.linspace(-0.5, 0.5, n_fine, endpoint=False)
    return xf, np.power(10.0, np.interp(xf, ph_ext, lb_ext))


# --------------------------------------------------------------------------- #
# Jitter (dual-Dirac) y convolución con la bañera
# --------------------------------------------------------------------------- #
def jitter_pdf(xf, sigma_rj, dj_pp, modelo="dual_dirac"):
    """PDF de jitter sobre la rejilla xf (UI, centrada en 0):
      - 'dual_dirac': RJ gaussiano (sigma) + DJ como dos Diracs en +/- dj_pp/2
        (cota de peor caso; corte abrupto).
      - 'pj': DJ periódico SENOIDAL -> distribución arcoseno sobre [-dj_pp/2, dj_pp/2]
        convolucionada con el RJ gaussiano (más realista; suaviza el acantilado)."""
    sg = max(sigma_rj, 1e-9)
    if modelo == "pj":
        A = dj_pp / 2.0
        arc = np.zeros_like(xf)
        if A > 1e-9:
            m = np.abs(xf) < A * 0.999
            arc[m] = 1.0 / (np.pi * np.sqrt(A ** 2 - xf[m] ** 2))
            arc = arc / (arc.sum() + 1e-18)
        else:
            arc[np.argmin(np.abs(xf))] = 1.0
        g = np.exp(-0.5 * (xf / sg) ** 2)
        g = g / (g.sum() + 1e-18)
        conv = np.fft.fftshift(np.real(np.fft.ifft(
            np.fft.fft(np.fft.ifftshift(arc)) * np.fft.fft(np.fft.ifftshift(g)))))
        p = np.clip(conv, 0.0, None)
    else:  # dual_dirac
        g = lambda mu: np.exp(-0.5 * ((xf - mu) / sg) ** 2)
        p = 0.5 * g(-dj_pp / 2.0) + 0.5 * g(+dj_pp / 2.0)
    s = p.sum()
    return p / s if s > 0 else p


def aplicar_jitter(xf, ber_fine, sigma_rj, dj_pp, modelo="dual_dirac"):
    """BER observada con jitter = convolución circular de la bañera con la PDF."""
    kern = jitter_pdf(xf, sigma_rj, dj_pp, modelo)
    K = np.fft.fft(np.fft.ifftshift(kern))     # centrar el kernel en 0
    B = np.fft.fft(ber_fine)
    out = np.real(np.fft.ifft(B * K))
    return np.clip(out, 1e-12, 1.0)


def eye_width(xf, ber, umbral=1e-2):
    """Apertura horizontal (UI): tramo contiguo en torno al mínimo con BER<=umbral."""
    nfb = len(ber)
    below = ber <= umbral
    if not below.any():
        return 0.0
    imin = int(np.argmin(ber))
    if not below[imin]:
        return 0.0
    l = imin
    while (imin - l) < nfb and below[(l - 1) % nfb]:
        l -= 1
    r = imin
    while (r - imin) < nfb and below[(r + 1) % nfb]:
        r += 1
    dx = xf[1] - xf[0]
    return float((r - l) * dx)


# --------------------------------------------------------------------------- #
# Jitter total dual-Dirac: TJ(BER) = DJ + 2 Q^{-1}(BER/2) sigma_RJ
# --------------------------------------------------------------------------- #
def tj_dual_dirac(ber, sigma_rj, dj_pp):
    """Jitter total pico-a-pico (UI) a una BER dada, modelo dual-Dirac."""
    from scipy.special import erfcinv
    q = np.sqrt(2.0) * erfcinv(2.0 * np.asarray(ber, dtype=float))   # Q^{-1}(BER)
    return dj_pp + 2.0 * q * sigma_rj


# --------------------------------------------------------------------------- #
# Envolvente VERTICAL del ojo (apertura del peor ojo en V) vs fase — contorno 2D
# --------------------------------------------------------------------------- #
def envolvente_vertical(y, sym, sps_val=SPS, base=None, k_sigma=3.0, n_descartar=100):
    """Devuelve (fase_ui, V) con la apertura del peor ojo PAM-4 (a k·sigma) en
    voltios, en función de la fase de muestreo dentro de la UI. Traza el contorno
    vertical del ojo 2D. (V está escalada por la ganancia del CTLE.)"""
    from cadena_enlace import estadisticas_niveles as _est
    if base is None:
        base = alinear_base(y, sym, sps_val)
    n = len(y)
    sym = np.asarray(sym, dtype=float)
    Nsym = len(sym)
    deltas = np.arange(-sps_val // 2, sps_val // 2)
    V = np.full(len(deltas), -99.0)
    for i, d in enumerate(deltas):
        off = base + d
        idx = off + sps_val * np.arange(Nsym)
        idx = idx[(idx >= 0) & (idx < n)]
        kk = len(idx)
        if kk <= n_descartar + 20:
            continue
        m = y[idx][n_descartar:]
        tx = sym[:kk][n_descartar:]
        mu, s, cnt = _est(m, tx)
        if np.any(cnt == 0):
            continue
        o = np.argsort(mu)
        mu, s = mu[o], s[o]
        ap = (mu[1:] - k_sigma * s[1:]) - (mu[:-1] + k_sigma * s[:-1])
        V[i] = float(np.min(ap))
    return deltas / sps_val, V

