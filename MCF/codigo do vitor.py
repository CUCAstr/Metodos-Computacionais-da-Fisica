import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as pla
import scipy.integrate as spi

# Parâmetros conforme solicitado
L0 = 3.0    # Distância entre pivôs e comprimento de repouso da mola (m)
L  = 2.0    # Comprimento dos pêndulos (m)
g  = 9.81   # Gravidade (m/s^2)
k  = 20.0   # Constante elástica (N/m)
m1 = 1.0    # Massa 1 (kg)
m2 = 1.0    # Massa 2 (kg)

def pendulos_ligados(Y, t):
    tht1, tht2, vt1, vt2 = Y
    
    # f = s^2 = (x2 - x1)^2 + (y2 - y1)^2
    # onde s é a distância entre as massas
    
    # Posições para o cálculo da distância da mola
    term_x = L0 + L * np.sin(tht2) - L * np.sin(tht1)
    term_y = L * np.cos(tht1) - L * np.cos(tht2)
    
    # f = s^2
    f = term_x**2 + term_y**2
    s = np.sqrt(f)
    eps = s - L0 # Deformação da mola
    
    # Derivadas de f em relação a theta1 e theta2
    df_dtht1 = -2 * L * term_x * np.cos(tht1) + 2 * L * term_y * (-np.sin(tht1))
    df_dtht2 =  2 * L * term_x * np.cos(tht2) + 2 * L * term_y * (np.sin(tht2))
    
    # Derivada da deformação: deps/dth = (1/2s) * df/dth
    deps_dtht1 = df_dtht1 / (2 * s)
    deps_dtht2 = df_dtht2 / (2 * s)
  
    dtht1 = vt1
    dtht2 = vt2
    
    # Equações de movimento (Torque / Momento de Inércia m*L^2)
    # dvt = [Torque_gravidade + Torque_mola] / (m*L^2)
    # Torque_mola = -d/dtheta (1/2 * k * eps^2) = -k * eps * deps/dtheta
    
    dvt1  = -(g/L) * np.sin(tht1) - (k * eps * deps_dtht1) / (m1 * L**2)
    dvt2  = -(g/L) * np.sin(tht2) - (k * eps * deps_dtht2) / (m2 * L**2)
    
    return [dtht1, dtht2, dvt1, dvt2]

# Condições iniciais
tht10, tht20, vt10, vt20 = np.pi/12, -np.pi/12, 0, 0
Y0 = [tht10, tht20, vt10, vt20]

# Tempo de simulação
tmin, tmax, h = 0, 20, 0.04
t = np.arange(tmin, tmax+h, h)

# Integração
Y = spi.odeint(pendulos_ligados, Y0, t)
tht1, tht2, vt1, vt2 = Y.T

# Coordenadas cartesianas para animação
x1 = L * np.sin(tht1)
y1 = -L * np.cos(tht1)
x2 = L0 + L * np.sin(tht2)
y2 = -L * np.cos(tht2)

# Configuração da animação
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-L, L0 + L), ylim=(-L - 1, 0.5))
ax.grid()

fio1, = ax.plot([], [], '-k', lw=2)
fio2, = ax.plot([], [], '-k', lw=2)
mola, = ax.plot([], [], '--r', lw=1.5)
massa1, = ax.plot([], [], 'ob', markersize=12)
massa2, = ax.plot([], [], 'og', markersize=12)

def update(n):
    fio1.set_data([0, x1[n]], [0, y1[n]])
    fio2.set_data([L0, x2[n]], [0, y2[n]])
    mola.set_data([x1[n], x2[n]], [y1[n], y2[n]])
    massa1.set_data([x1[n]], [y1[n]])
    massa2.set_data([x2[n]], [y2[n]])
    return fio1, fio2, mola, massa1, massa2

ani = pla.FuncAnimation(fig, update, frames=len(t), interval=40, blit=True)
plt.show()
