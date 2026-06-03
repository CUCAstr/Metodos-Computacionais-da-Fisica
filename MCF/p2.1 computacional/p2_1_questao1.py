import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt

'''
LOM3227 — Métodos Computacionais da Física
Resolução da P2.1c - Questão 1
'''

def f_temporal(t):
    """Função f(t) conforme definida na Eq (1) da prova"""
    res = np.zeros_like(t)
    res[(t >= -1) & (t < 0)] = -1
    res[(t >= 0) & (t <= 1)] = 1
    return res

def F_exata(w):
    """Transformada de Fourier analítica fornecida"""
    # Evitar divisão por zero em w=0 (limite é 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        res = -1j * np.sqrt(2/np.pi) * (1 - np.cos(w)) / w
    res[w == 0] = 0
    return res

# Parâmetros para a FFT (Dica 1: período de observação bem grande)
T_obs = 60.0        
N_fft = 2**15       # Alta amostragem para minimizar aliasing e dispersão
dt = T_obs / N_fft

# Vetor de tempo centrado para amostrar a função f(t) definida em [-1, 1]
t_fft = np.linspace(-T_obs/2, T_obs/2, N_fft, endpoint=False)

# Amostragem do sinal
sinal = f_temporal(t_fft)

# Cálculo da FFT
# Rotacionamos o sinal para que t=0 seja a primeira posição (convenção FFT)
sinal_shifted = fft.ifftshift(sinal)
FD_w = fft.fft(sinal_shifted)
freqs = fft.fftfreq(N_fft, d=dt)
w_discreta = 2 * np.pi * freqs

# Relação teórica vs discreta (Dica 3)
# F(w) = (N/T) * (1/sqrt(2pi)) * FD(w)
F_discreta = (N_fft / T_obs) * (1 / np.sqrt(2 * np.pi)) * FD_w

# Filtragem para o intervalo solicitado: 0 <= w <= 6pi
mask = (w_discreta >= 0) & (w_discreta <= 6 * np.pi)
w_plot = w_discreta[mask]
abs_F_disc = np.absolute(F_discreta[mask])

# Cálculo da solução exata para comparação
w_cont = np.linspace(0.001, 6 * np.pi, 1000)
abs_F_exata = np.absolute(F_exata(w_cont))

# Normalização (Dica 2: máximos parecidos)
abs_F_disc /= np.max(abs_F_disc)
abs_F_exata /= np.max(abs_F_exata)

# Gráfico
plt.figure(figsize=(10, 6))
plt.plot(w_cont, abs_F_exata, 'r-', lw=2, label='Analítica $|F(\omega)|$')
plt.plot(w_plot, abs_F_disc, 'b--', label='FFT $|F(\omega)|$ (Normalizada)')
plt.title('Questão 1: Transformada de Fourier Exata vs Discreta')
plt.xlabel('$\omega$ (rad/s)')
plt.ylabel('Módulo Normalizado')
plt.xlim(0, 6 * np.pi)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
