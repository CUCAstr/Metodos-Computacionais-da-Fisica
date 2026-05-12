import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

L = 8  # m
rho0 = 4  # kg/m
g = 9.81  # m/s²
T = 5000  # N

c2 = T / rho0  # (m/s)²

# discretização da posição
Nx = 250
x = np.linspace(0, L, Nx+1)
Dx = L / Nx

# discretização do tempo
Nt = 2 ** 17
tmax = 100  # s
Dt = tmax / Nt
t = np.linspace(0, tmax, Nt+1)

# declaração do array deslocamento
u = np.zeros((Nt+1, Nx+1))

# condições iniciais
u0 = np.zeros(Nx+1)  # formato inicial da corda (horizontal)
# u0 = .05 * np.sin(3 * x *np.pi / L)
v0 = np.zeros(Nx+1)  # velocidade inicial da corda (repouso)

# Diferenças finitas

u[0, :] = np.copy(u0)  # formato inicial da corda

f = c2 * (Dt / Dx)**2

# Deslocamento no instante t1:
u[1, 1:-1] = u[0, 1:-1] + Dt * v0[1:-1] + .5 * f * (u[0, 2:] - 2 * u[0, 1:-1] + u[0, :-2]) - .5 * g * Dt**2

# Deslocamento nos próximos instantes tj:
for j in range(1, Nt):
    u[j+1, 1:-1] = f * (u[j, :-2] + u[j, 2:]) + 2 * (1 - f) * u[j, 1:-1] - u[j-1, 1:-1] - g * Dt ** 2

fig = plt.figure()

plt.xlim(0, L)
plt.ylim(u.min(), u.max())

plt.plot(x, u0, '--')

corda, = plt.plot([], [])

def update(j):
    corda.set_data(x, u[j, :])
    return corda,

a = anim.FuncAnimation(fig, update, blit=True, frames=Nt, interval=5)

plt.show()