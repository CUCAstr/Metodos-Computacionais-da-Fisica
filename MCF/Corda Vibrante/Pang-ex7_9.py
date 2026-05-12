import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import matplotlib.pyplot as plt
import numpy.fft as fft
import matplotlib.animation as anim

# Dados (Pang, Ex. 7.9)

L = 8.0
T = 5.0e3
g = 9.81

Nx = 250
x = np.linspace(0, L, Nx + 1)
h = L / Nx

def amigo(mi, xi, ai, x):
    gauss = np.zeros_like(x)
    trunc = np.where(np.abs(x - xi) <= ai)
    gauss[trunc] = .5 * mi / ai * np.exp(-((x[trunc] - xi) / ai)**2)
    return gauss

corda = 4 * np.ones_like(x)
amigo1 = amigo(40, 3.5, 0.3, x)
amigo2 = amigo(55, 4.1, 0.27, x)

densidade = corda + amigo1 + amigo2
c = np.sqrt(T / densidade)

# equilíbrio

P = g * (h / c) ** 2
MI = sps.diags([np.ones(Nx-2), -2 * np.ones(Nx-1), np.ones(Nx-2)], [-1, 0, 1]).toarray()
U0 = np.zeros_like(x)
U0[1:-1] = la.solve(MI, P[1:-1])  # deslocamento de equilíbrio; as extremidades ficam em repouso

# Frequências de vibração

M = sps.diags([-np.ones(Nx-2), 2 * np.ones(Nx-1), -np.ones(Nx-2)], [-1, 0, 1]).toarray()
R = np.diag((h / c[1:-1])**2)
autoval, autovec = la.eigh(M, R)
freq = np.sqrt(autoval)  # Frequências de vibração
print('primeiras cinco frequências (Hz):', freq[:5])

# Solução usando diferenças finitas

u0 = np.zeros_like(x)  # condição inicial
# u0 = .05 * np.sin(3 * x *np.pi / L)
du0 = np.zeros_like(x)  # derivada espacial na condição inicial

tmin, tmax, Nt = 0, 100, 2**17
t = np.linspace(tmin, tmax, Nt+1)
dt = t[1] - t[0]
u = np.zeros((t.size, x.size))
f = (dt * c / h) ** 2

u[0, :] = np.copy(u0)
u[1, 1:-1] = .5 * f[1:-1] * (u[0, :-2] + u[0, 2:]) + (1 - f[1:-1]) * u[0, 1:-1] + dt * du0[1:-1] - .5 * g * dt ** 2

for i in range(2, t.size):
    u[i, 1:-1] = f[1:-1] * (u[i-1, :-2] + u[i-1, 2:]) + 2 * (1 - f[1:-1]) * u[i-1, 1:-1] - u[i-2, 1:-1] - g * dt ** 2

# Frequências do ponto central da corda - via FFT

wNyq = t.size * np.pi / tmax
dw = 2 * wNyq / t.size
w = np.linspace(0, 2 * wNyq, t.size)
print(f'w_Nyq = {wNyq:.3f} Hz, dw = {dw:.3f} Hz')

fft = fft.fft(u[:, x.size // 2])
power = np.absolute(fft) ** 2
power /= power.max() / 100  # normalizando em [0-100]

# Gráficos e animação

# Animação do movimento

fig = plt.figure(figsize=(12, 6))
ax = plt.axes()
ax2 = ax.twinx()

ax.set_xlabel('$x$ (m)')
ax.set_ylabel('deslocamento (m)')
ax2.set_ylabel('densidade (kg/m)')

ax.plot(x, U0, 'r--', label='equilíbrio')
ax.plot(x, u0, 'r-.', label='cond. inicial')
ax2.plot(x, densidade, 'b', label='densidade')

corda, = ax.plot(x, u0, 'r', label='$u(x, t)$')
tempo = ax.text(.1, .1, '', transform=ax.transAxes)

ax.set_xlim(0, L)
ax.set_ylim(u.min(), u.max())


def update(i):
    corda.set_xdata(x)
    corda.set_ydata(u[i, :])
    tempo.set_text(f'$t = {t[i]:.2f}$ s')
    return corda, tempo


a = anim.FuncAnimation(fig, update, frames=t.size, interval=1, blit=True)
ax.legend(loc=2)
ax2.legend()
plt.tight_layout()

# gráfico da FFT

# plt.figure()
#
# plt.title('Transformada de Fourier do deslocamento central')
# plt.xlabel('$\\omega$ (Hz)')
# plt.ylabel('Potência espectral (u.a.)')
# plt.yscale('log')
#
# wmax = freq[5]
# plt.plot(w[w < wmax], power[w < wmax])
#
# for i in range(5):
#     plt.axvline(freq[i], color='r', alpha=.3)

plt.show()
