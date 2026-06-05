import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.fft import fft, fftfreq

# ==========================================
# Parâmetros do Problema
# ==========================================
a = 1.3          # m
w0 = 1.5         # s^-1 (rad/s)
g = 9.81         # m/s^2
L = 0.7          # m
m = 0.5          # kg (não é usado na EDO final, mas listado por completude)
theta0 = np.pi/6 # rad
dtheta0 = 0.0    # "não há movimento da massa em relação ao carrinho" no t=0

# ==========================================
# Parte (a): EDO e Integração Numérica
# ==========================================
def pendulum_ode(t, y):
    theta, omega = y
    # y[0] = theta, y[1] = dtheta/dt (omega)
    dydt = [
        omega,
        -(g / L) * np.sin(theta) + (a * w0**2 / L) * np.cos(w0 * t) * np.cos(theta)
    ]
    return dydt

# Intervalo para o item (a)
t_span_a = (0, 10)
t_eval_a = np.linspace(0, 10, 1000)

sol_a = solve_ivp(pendulum_ode, t_span_a, [theta0, dtheta0], t_eval=t_eval_a, method='RK45')

plt.figure(figsize=(10, 4))
plt.plot(sol_a.t, sol_a.y[0], 'b-', linewidth=2)
plt.title(r'Item (a): Evolução temporal de $\theta(t)$ para $0 \leq t \leq 10$ s')
plt.xlabel('Tempo (s)')
plt.ylabel(r'Ângulo $\theta$ (rad)')
plt.grid(True)
plt.tight_layout()
plt.show()

# ==========================================
# Parte (b): Transformada Discreta de Fourier
# ==========================================
T_obs = 200.0
N_pts = 2**13  # 8192 pontos
t_eval_b = np.linspace(0, T_obs, N_pts, endpoint=False)
dt = T_obs / N_pts

# Integração para 200s garantindo a amostragem exata
sol_b = solve_ivp(pendulum_ode, (0, T_obs), [theta0, dtheta0], t_eval=t_eval_b, method='RK45')
theta_t = sol_b.y[0]

# Cálculo da FFT
theta_fft = fft(theta_t)
freqs = fftfreq(N_pts, dt)
w = 2 * np.pi * freqs  # Convertendo para frequência angular (rad/s)

# Frequências esperadas para comparação
wp = np.sqrt(g / L)  # Frequência natural do pêndulo

# Pegando apenas as frequências positivas para o gráfico
positive_idx = w > 0
w_pos = w[positive_idx]
mag_fft = np.abs(theta_fft[positive_idx]) / (N_pts / 2) # Normalizando a magnitude

plt.figure(figsize=(10, 5))
plt.plot(w_pos, mag_fft, 'k-', linewidth=1.5)

# Marcando as frequências w0 e wp
plt.axvline(w0, color='r', linestyle='--', label=r'Frequência de Forçamento ($\omega_0 = 1.50$ rad/s)')
plt.axvline(wp, color='g', linestyle='--', label=rf'Frequência Natural ($\omega_p \approx {wp:.2f}$ rad/s)')

plt.xlim(0, 8) # Limitando o eixo x para focar nas frequências de interesse
plt.title(r'Item (b): Espectro de Frequências (FFT) de $\theta(t)$')
plt.xlabel(r'Frequência Angular $\omega$ (rad/s)')
plt.ylabel(r'Magnitude da Transformada $|\Theta(\omega)|$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()