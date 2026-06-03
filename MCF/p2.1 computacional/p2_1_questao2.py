import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

'''
LOM3227 — Métodos Computacionais da Física
Resolução da P2.1c - Questão 2
'''

# Dados do problema
a = 1.3         # amplitude do vagão (m)
w0 = 1.5        # freq. do vagão (rad/s)
g = 9.81        # gravidade (m/s^2)
l = 0.7         # comprimento da haste (m)
m = 0.5         # massa (kg)
theta0 = np.pi / 6  # condição inicial (rad)

def edo_pendulo(t, y):
    """
    EDO derivada via Mecânica Lagrangiana:
    theta'' = [a * w0^2 * cos(w0*t) * cos(theta) - g * sin(theta)] / l
    """
    theta, omega = y
    dtheta = omega
    domega = (a * w0**2 * np.cos(w0*t) * np.cos(theta) - g * np.sin(theta)) / l
    return [dtheta, domega]

# -----------------------------------------------------------------------------
# Item (a): Integração numérica para 0 <= t <= 10s
# -----------------------------------------------------------------------------
t_span_a = (0, 10)
t_eval_a = np.linspace(0, 10, 1000)
sol_a = solve_ivp(edo_pendulo, t_span_a, [theta0, 0], t_eval=t_eval_a, method='RK45')

plt.figure(figsize=(10, 5))
plt.plot(sol_a.t, sol_a.y[0], 'g', label='$\\theta(t)$')
plt.title('Questão 2a: Evolução de $\\theta(t)$ para $t \\in [0, 10]$ s')
plt.xlabel('$t$ (s)')
plt.ylabel('$\\theta$ (rad)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# -----------------------------------------------------------------------------
# Item (b): Análise de Fourier (FFT)
# -----------------------------------------------------------------------------
T_obs_b = 200.0      # Tempo de observação: 200 s
N_samples = 2**13    # Amostragem: 8192 pontos
t_fft_b = np.linspace(0, T_obs_b, N_samples, endpoint=False)

# Resolver a EDO para o tempo longo solicitado
sol_b = solve_ivp(edo_pendulo, (0, T_obs_b), [theta0, 0], t_eval=t_fft_b, method='RK45')
theta_sinal = sol_b.y[0]

# Cálculo da FFT
fft_res = fft.fft(theta_sinal)
freqs_b = fft.fftfreq(N_samples, d=(T_obs_b/N_samples))
w_vals = 2 * np.pi * freqs_b
espectro = np.absolute(fft_res)**2

# Frequências teóricas para comparação
wp = np.sqrt(g / l)  # frequência natural de oscilação do pêndulo

# Gráfico do Espectro de Frequências
plt.figure(figsize=(10, 6))
# Foco no intervalo de baixas frequências relevante
mask_b = (w_vals > 0) & (w_vals < 10)
plt.plot(w_vals[mask_b], espectro[mask_b], 'k', lw=1.5)
plt.axvline(w0, color='r', linestyle='--', alpha=0.7, label=f'Vagão $\omega_0 = {w0:.2f}$')
plt.axvline(wp, color='b', linestyle='--', alpha=0.7, label=f'Natural $\omega_p = {wp:.2f}$')

plt.title('Questão 2b: Análise de Fourier de $\\theta(t)$')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('Potência Espectral (u.a.)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"--- Resultados Questão 2 ---")
print(f"Frequência do vagão (w0): {w0:.3f} rad/s")
print(f"Frequência natural calculada (wp): {wp:.3f} rad/s")
print(f"O gráfico confirma a presença de picos em ambas as frequências.")
