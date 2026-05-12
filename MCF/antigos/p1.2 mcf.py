import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.animation as anim

g = 9.81  # m/s²
M = 2.0   # kg (massa do carrinho)
m = 1.0   # kg (massa do pêndulo)
K = 50.0  # N/m (mola da parede)
k = 40.0  # N/m (mola do pêndulo)
L0 = 1.0  # m (comprimento natural)

# sistema de EDOS de 1ª ordem:
def carrinho_pendulo(Y, t):
    eta, eps, theta, veta, veps, vtheta = Y
    
    deta = veta
    deps = veps
    dtheta = vtheta
    
    # Matriz M e Vetor B para resolver as acelerações (dveta, dveps, dvtheta)
    M_mat = np.array([
        [M + m, m * np.sin(theta), m * (L0 + eps) * np.cos(theta)],
        [m * np.sin(theta), m, 0],
        [m * (L0 + eps) * np.cos(theta), 0, m * (L0 + eps)**2]
    ])
    
    B1 = -K * eta - 2 * m * veps * vtheta * np.cos(theta) + m * (L0 + eps) * (vtheta**2) * np.sin(theta)
    B2 = -k * eps + m * (L0 + eps) * (vtheta**2) + m * g * np.cos(theta)
    B3 = -2 * m * (L0 + eps) * veps * vtheta - m * g * (L0 + eps) * np.sin(theta)
    B = np.array([B1, B2, B3])
    
    # Resolve o sistema linear M * a = B
    dveta, dveps, dvtheta = np.linalg.solve(M_mat, B)
    
    return [deta, deps, dtheta, dveta, dveps, dvtheta]

# condição inicial:
eta0, eps0, theta0 = 0.5, 0.0, 0.3    # m, m, rad
veta0, veps0, vtheta0 = 0.0, 0.0, 0.0 # m/s, m/s, rad/s
Y0 = [eta0, eps0, theta0, veta0, veps0, vtheta0]

# tempo:
tmin, tmax, h = 0, 20, 0.033
t = np.arange(tmin, tmax+h, h)

# solução com odeint:
Y = spi.odeint(carrinho_pendulo, Y0, t)

# separando as incógnitas:
eta, eps, theta, veta, veps, vtheta = Y.T

# posição do carrinho e da massa m:
xp = eta
yp = np.zeros_like(eta)

xm = eta + (L0 + eps) * np.sin(theta)
ym = - (L0 + eps) * np.cos(theta)

# gráfico:
fig = plt.figure()

# Limites dinâmicos
plt.ylim(ym.min() - 0.5, 1)
plt.xlim(min(xp.min(), xm.min()) - 1, max(xp.max(), xm.max()) + 1)
plt.grid(True, linestyle='--', alpha=0.6)

# Inicialização dos elementos
mola_K, = plt.plot([], [], '--k', label='Mola K (Parede)')       # mola da parede
mola_k, = plt.plot([], [], '-r', label='Mola k (Pêndulo)')       # mola do pêndulo
carrinho, = plt.plot([], [], 'sk', ms=8, label='Carrinho M')     # quadrado preto (carrinho)
bola, = plt.plot([], [], 'ob', ms=10, label='Massa m')           # bola azul (massa m)
traj, = plt.plot([], [], '-b', alpha=0.3, label='Trajetória')    # trajetória 

plt.legend(loc='upper right', fontsize=9)

# função de atualização:
def update(n):
    # Parede fixa num x arbitrário à esquerda (x = -3)
    mola_K.set_data([-3, xp[n]], [0, 0]) 
    mola_k.set_data([xp[n], xm[n]], [yp[n], ym[n]])
    carrinho.set_data([xp[n]], [yp[n]])
    bola.set_data([xm[n]], [ym[n]])
    traj.set_data(xm[:n], ym[:n])
    
    return mola_K, mola_k, carrinho, bola, traj,

# controle da animação:
a = anim.FuncAnimation(fig, update, blit=True, frames=t.size, interval=h*1000)

plt.show()