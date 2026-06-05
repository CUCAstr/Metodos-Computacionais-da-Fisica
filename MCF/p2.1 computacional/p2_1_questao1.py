import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt

def f(t):
    # f(t) da questao 1
    y = np.zeros_like(t)
    y[(t >= -1) & (t < 0)] = -1
    y[(t >= 0) & (t <= 1)] = 1
    return y

def F_analitica(w):
    # formula da transformada exata dada no pdf
    res = np.zeros_like(w, dtype=complex)
    m = (w != 0)
    res[m] = -1j * np.sqrt(2/np.pi) * (1 - np.cos(w[m])) / w[m]
    return res

# parametros
T = 100
N = 2**16
dt = T/N
t = np.linspace(-T/2, T/2, N, endpoint=False)

sinal = f(t)

# fft
sinal_p = fft.ifftshift(sinal)
FD = fft.fft(sinal_p)
freqs = fft.fftfreq(N, d=dt)
w_d = 2 * np.pi * freqs

# correcao de escala
F_d = (N / (T * np.sqrt(2 * np.pi))) * FD

# plot
mask = (w_d >= 0) & (w_d <= 6 * np.pi)
w_p = w_d[mask]
y_d = np.absolute(F_d[mask])

w_c = np.linspace(0.01, 6 * np.pi, 500)
y_a = np.absolute(F_analitica(w_c))

# normalizando
y_d = y_d / np.max(y_d)
y_a = y_a / np.max(y_a)

plt.figure()
plt.plot(w_c, y_a, 'r-', label='Exata')
plt.plot(w_p, y_d, 'b.', alpha=0.5, label='FFT')
plt.title('Questao 1')
plt.xlabel('omega')
plt.ylabel('|F(omega)|')
plt.legend()
plt.grid()
plt.show()
