# %%
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.animation as anim


# %%
# Parâmetros do sistema
g = 9.81   # m/s^2
m1 = 1.0   # kg
m2 = 1.0   # kg
k = 20.0   # N/m
L = 2.0    # Comprimento dos pêndulos (m)
L0 = 3.0   # Distância entre pivôs e comp. de repouso da mola (m)

# Sistema de EDOs de 1ª ordem:
def pendulos_acoplados(Y, t):
    th1, th2, v1, v2 = Y
    
    # Diferenças de posição geométrica
    dx = L0 + L * (np.sin(th2) - np.sin(th1))
    dy = L * (np.cos(th2) - np.cos(th1))
    
    # Distância atual entre as massas
    s = np.sqrt(dx**2 + dy**2)
    
    # Prevenção de divisão por zero na mola
    if s == 0:
        fm = 0
    else:
        fm = k * (s - L0) / s
        
    # Derivadas
    dth1 = v1
    dth2 = v2
    dv1 = -(g/L)*np.sin(th1) + (fm/(m1*L)) * (dx*np.cos(th1) - dy*np.sin(th1))
    dv2 = -(g/L)*np.sin(th2) - (fm/(m2*L)) * (dx*np.cos(th2) - dy*np.sin(th2))
    
    return [dth1, dth2, dv1, dv2]

# Condições iniciais (th1, th2, v1, v2)
# Exemplo: Pêndulo 1 inclinado, Pêndulo 2 em repouso
Y0 = [0.5, 0.0, 0.0, 0.0] 

# Tempo de simulação:
tmin, tmax, h = 0, 20, 0.04
t = np.arange(tmin, tmax+h, h)

# Solução com odeint:
Y = spi.odeint(pendulos_acoplados, Y0, t)

# Separando as incógnitas:
th1, th2, v1, v2 = Y.T

# Posições cartesianas para a animação (invertendo Y para visualização correta)
x1 = L * np.sin(th1)
y1 = -L * np.cos(th1)

x2 = L0 + L * np.sin(th2)
y2 = -L * np.cos(th2)

# Gráfico e Animação:
fig = plt.figure(figsize=(8, 5))
plt.xlim(-L, L0 + L)
plt.ylim(-L - 1, 0.5)
plt.grid(True)
plt.gca().set_aspect('equal')

# Inicialização das linhas e pontos
fio1, = plt.plot([], [], '-k', lw=2)
fio2, = plt.plot([], [], '-k', lw=2)
mola, = plt.plot([], [], '--r', lw=1.5)
massa1, = plt.plot([], [], 'ob', markersize=10)
massa2, = plt.plot([], [], 'og', markersize=10)

# Função de atualização
def update(n):
    fio1.set_data([0, x1[n]], [0, y1[n]])
    fio2.set_data([L0, x2[n]], [0, y2[n]])
    mola.set_data([x1[n], x2[n]], [y1[n], y2[n]])
    massa1.set_data([x1[n]], [y1[n]])
    massa2.set_data([x2[n]], [y2[n]])
    return fio1, fio2, mola, massa1, massa2

# Controle da animação:
a = anim.FuncAnimation(fig, update, frames=t.size, interval=40, blit=True)

plt.show()

