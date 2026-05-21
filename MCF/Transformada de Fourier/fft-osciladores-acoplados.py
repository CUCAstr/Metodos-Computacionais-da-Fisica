import numpy as np
import scipy.integrate as spi
import scipy.fftpack as fourier
import matplotlib.pyplot as plt

k, K = 1, 1.5
m = 0.3

x10, x20, v10, v20 = 0.2, 0.5, 0.8, 0
Y0 = [x10, x20, v10, v20]

w1q = k / m
w2q = K / m

# Sistemas de EDOs
def osc_harm_acop(Y, t):
	x1, x2, v1, v2 = Y
	dx1 = v1
	dx2 = v2
	dv1 = -(w1q + w2q) * x1 + w2q * x2
	dv2 = w2q * x1 - (w1q + w2q) * x2
	return [dx1, dx2, dv1, dv2]

T = 200 # periodo de observacao
N = 2 ** 12 # amostragem

t = np.linspace(0, T, N+1)

# Solucao das EDOs
Y = spi.odeint(osc_harm_acop, Y0, t)
x1, x2, v1, v2 = Y.T

#  Transformada de fourier de x1
F = fourier.fft(x1)
Fa = np.absolute(F)

wNyq = np.pi * N / T
w = np.linspace(0, 2 * wNyq, N + 1)

# frequencias teoricas
wa = np.sqrt(w1q)
wb = np.sqrt(w1q + 2 * w2q)

# Graficos
plt.figure()
plt.plot(w, Fa)

plt.axvline(wa, color='r')
plt.axvline(wb, color='b')
plt.xlim(0, wNyq)

plt.show()


