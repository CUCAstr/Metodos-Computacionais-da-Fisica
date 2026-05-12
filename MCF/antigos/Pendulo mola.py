import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.animation as anim

g = 9.81  # m/s²
m = 0.3  # kg
k = 50  # N / m
L0 = 3  # m

# sistema de EDOS de 1ª ordem:
def pendulo_elastico(Y, t):
    eps, theta, veps, vtheta = Y
    deps = veps
    dtheta = vtheta
    dveps = (L0 + eps) * vtheta ** 2 + g * np.cos(theta) - k * eps / m
    dvtheta = -(g * np.sin(theta) + 2 * veps * vtheta) / (L0 + eps)
    return [deps, dtheta, dveps, dvtheta]

# condição inicial:
eps0, theta0, veps0, vtheta0 = 0.3, np.pi/3, 0, 0  # m, rad, m/s, rad/s
Y0 = [eps0, theta0, veps0, vtheta0]

# tempo:
tmin, tmax, h = 0, 20, 0.01
t = np.arange(tmin, tmax+h, h)

# solução com odeint:
Y = spi.odeint(pendulo_elastico, Y0, t)

# separando as incógnitas:
eps, theta, veps, vtheta = Y.T

# posição da massa m:
x = (L0 + eps) * np.sin(theta)
y = - (L0 + eps) * np.cos(theta)

# gráfico:
fig = plt.figure()

ymax = 0 if y.max() < 0 else y.max()
plt.ylim(y.min(), ymax)
plt.xlim(x.min(), x.max())

bola, = plt.plot([], [], 'o-r')  # inicialização da bola
traj, = plt.plot([], [], '-b')  # inicialização da trajetória

# função de atualização:
def update(n):
    bola.set_data([0, x[n]], [0, y[n]])
    traj.set_data(x[:n], y[:n])
    return bola, traj,

# controle da animação:
a = anim.FuncAnimation(fig, update, blit=True, frames=t.size, interval=15)

plt.show()