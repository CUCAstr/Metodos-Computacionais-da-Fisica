import numpy as np
import scipy.integrate as spi
import scipy.fftpack as fourier
import matplotlib.pyplot as plt

# Sistemas de EDOs
def vanderpol(Y, t):
	x, v = Y
	dx = v
	dv = -x - eps * (x**2 - 1) * v
	return [dx, dv]

eps = 1
tmax = 8 * np.pi
N = 2**10
t = np.linspace(0, tmax, N+1)
T = 200

x0, v0 = 0.5, 0
Y0 = [x0, v0]

Y = spi.odeint(vanderpol, Y0, t)
x, v = Y.T

# Transformada de fourier de x1
F = fourier.fft(x)
Fa = np.absolute(F)

wNyq = np.pi * N / T
w = np.linspace(0, 2 * wNyq, N + 1)


# Graficos
plt.figure()
plt.plot(w, Fa)
plt.xlim(0, wNyq)

plt.show()


