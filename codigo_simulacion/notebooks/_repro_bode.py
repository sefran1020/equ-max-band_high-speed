import matplotlib
matplotlib.use('Agg')

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# Configuración estética de las gráficas
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['font.size'] = 10

# Parámetros típicos de una línea RC distribuida en una PCB de alta velocidad
# Longitud d = 15 cm
R_total = 100.0   # Resistencia total de la línea (Ohms)
C_total = 30e-12  # Capacitancia total de la línea (30 pF)

print(f"Parámetros de la línea:")
print(f"  R_total = {R_total} Ohms")
print(f"  C_total = {C_total * 1e12} pF")
print(f"  Constante de tiempo de difusión teórica Rt*Ct = {R_total * C_total * 1e9:.3f} ns")

def obtener_espacio_estados_rc(R_t, C_t, N):
    """
    Construye las matrices del espacio de estados (A, B, C, D)
    para una línea RC distribuida de N celdas en escalera.
    """
    R_seg = R_t / N
    C_seg = C_t / N
    alpha = 1.0 / (R_seg * C_seg)
    
    A = np.zeros((N, N))
    for i in range(N-1):
        A[i, i] = -2 * alpha
        A[i, i+1] = alpha
        A[i+1, i] = alpha
    # Condición de circuito abierto en el extremo de la línea
    A[N-1, N-1] = -1 * alpha
    
    B = np.zeros((N, 1))
    B[0, 0] = alpha
    
    C = np.zeros((1, N))
    C[0, N-1] = 1.0
    
    D = np.zeros((1, 1))
    
    return A, B, C, D

def evaluar_frecuencia_ss(A, B, C, D, w):
    """
    Evalúa la respuesta en frecuencia H(jw) = C * (jw*I - A)^-1 * B + D
    de forma directa y numéricamente estable.
    """
    N = A.shape[0]
    I = np.eye(N)
    H = []
    for omega in w:
        s = 1j * omega
        try:
            # Resolvemos (s*I - A) * x = B para hallar x = (s*I - A)^-1 * B
            x = np.linalg.solve(s * I - A, B)
            h = np.dot(C, x)[0, 0] + D[0, 0]
            H.append(h)
        except np.linalg.LinAlgError:
            H.append(np.nan)
    H = np.array(H)
    mag_db = 20 * np.log10(np.abs(H))
    fase_deg = np.rad2deg(np.unwrap(np.angle(H)))
    return mag_db, fase_deg

# Definimos el rango de frecuencias de análisis (10 MHz a 50 GHz)
freqs = np.logspace(7, 10.7, 1000)
w = 2 * np.pi * freqs

# 1. Solución analítica exacta de la línea de parámetros distribuidos
s = 1j * w
H_analitico = 1.0 / np.cosh(np.sqrt(s * R_total * C_total))
fase_analitico_deg = np.rad2deg(np.unwrap(np.angle(H_analitico)))

# Graficar
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

# Graficamos la solución analítica como línea negra gruesa de referencia
ax1.plot(freqs, 20 * np.log10(np.abs(H_analitico)), label='Analytic (Theoretical)', color='black', linewidth=3, zorder=10)
ax2.plot(freqs, fase_analitico_deg, color='black', linewidth=3, zorder=10)

# Graficamos aproximaciones para diferentes órdenes N
ordenes = [1, 2, 5, 15, 30]
colores = ['#e53e3e', '#dd6b20', '#4a5568', '#3182ce', '#319795']

for N, col in zip(ordenes, colores):
    A, B, C, D = obtener_espacio_estados_rc(R_total, C_total, N)
    mag, fase = evaluar_frecuencia_ss(A, B, C, D, w)
    
    ax1.plot(freqs, mag, label=f'Ladder N={N}', color=col, linestyle='--', linewidth=1.5)
    ax2.plot(freqs, fase, color=col, linestyle='--', linewidth=1.5)

ax1.set_xscale('log')
ax1.set_ylabel('Magnitude (dB)', fontsize=10)
ax1.legend()
ax1.set_ylim([-60, 5])

ax2.set_xscale('log')
ax2.set_ylabel('Phase (degrees)', fontsize=10)
ax2.set_xlabel('Frequency (Hz)', fontsize=10)
ax2.set_ylim([-270, 10])

plt.tight_layout()
fig.savefig('fig_bode.png', dpi=200, bbox_inches='tight')
print('SAVED fig_bode.png')
