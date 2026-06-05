import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift

# ==========================================
# 1. Parâmetros de Amostragem
# ==========================================
T = 100.0        # Período de observação muito grande para capturar a função inteira (de -1 a 1)
N = 2**16        # Número elevado de pontos (potência de 2 é ideal para FFT)
dt = T / N       # Passo de tempo (resolução temporal)

# Vetor de tempo de -T/2 a T/2
t = np.linspace(-T/2, T/2, N, endpoint=False)

# ==========================================
# 2. Definição da Função f(t)
# ==========================================
f_t = np.zeros_like(t)
f_t[(t >= -1) & (t < 0)] = -1
f_t[(t >= 0)  & (t <= 1)] = 1

# ==========================================
# 3. Transformada Discreta de Fourier (FFT)
# ==========================================
F_D = fft(f_t)

# Vetor de frequências cíclicas (f) e angulares (w = 2*pi*f)
freqs = fftfreq(N, dt)
w_D = 2 * np.pi * freqs

# Reorganizando os arrays para que a frequência 0 fique no centro
w_D_shifted = fftshift(w_D)
F_D_shifted = fftshift(F_D)

# Aplicando a relação teórica para converter a FFT discreta na contínua (Dica 3)
# F(w) = (T / N) * (1 / sqrt(2*pi)) * F_D(w)
F_FFT = (T / N) * (1 / np.sqrt(2 * np.pi)) * F_D_shifted
mag_F_FFT = np.abs(F_FFT)

# ==========================================
# 4. Transformada Analítica Exata
# ==========================================
def exata_F(w):
    # Evita divisão por zero substituindo w=0 por um valor ínfimo temporário
    w_safe = np.where(w == 0, 1e-10, w)
    return -1j * np.sqrt(2 / np.pi) * ((1 - np.cos(w_safe)) / w_safe)

mag_F_exata = np.abs(exata_F(w_D_shifted))
# Em w=0 a função é 0, o artifício do 1e-10 lida com isso perfeitamente

# ==========================================
# 5. Visualização (Gráfico)
# ==========================================
plt.figure(figsize=(10, 6))

# Plotando de 0 a 6*pi
plt.plot(w_D_shifted, mag_F_exata, 'r-', linewidth=3, label='Exata $|F(\omega)|$')
plt.plot(w_D_shifted, mag_F_FFT, 'b--', linewidth=2, label='FFT $|F_{FFT}(\omega)|$')

# Configurações do eixo X (0 <= w <= 6pi)
plt.xlim(0, 6 * np.pi)

# Para deixar o eixo X com ticks em múltiplos de pi
pi_ticks = [0, 2*np.pi, 4*np.pi, 6*np.pi]
pi_labels = ['0', r'$2\pi$', r'$4\pi$', r'$6\pi$']
plt.xticks(pi_ticks, pi_labels)

plt.title('Comparação: Transformada de Fourier Exata vs Discreta (FFT)')
plt.xlabel(r'Frequência Angular $\omega$ (rad/s)')
plt.ylabel(r'Módulo da Transformada $|F(\omega)|$')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()
plt.show()