import numpy as np
import scipy.linalg as la
import scipy.integrate as spi
import matplotlib.pyplot as plt

k = 10 # N/m
g = 9.81 # m/s^2
m = 2 # kg
L = 3 # m
L0 = 1 # m

def EDOs(Y,t):
    x, th, vx, vth = Y
    
    # encontrando dvx e dvth das Eqs de Euler-Lagrange:
    M = np.array([[1,         -L * np.cos(th)],
                  [np.cos(th), L             ]])
    
    b = np.array([L * vth**2 * np.sin(th) - 2 * k / m * x,
                  -g * np.sin(th)							]) 

    dvx, dvth = la.solve(M, b)
    
    # derivadas definidas pelas variaveis auxiliares (vx e vth) 
    dx = vx
    dth = vth
     
    return [dx, dth, dvx, dvth]

# Condicoes iniciais

x0, th0, vx0, vth0 = 1, np.radians(30), 0, 0
Y0 = [x0, th0, vx0, vth0]

# Tempo
tmax, dt = 60, 0.01
t = np.arange(0, tmax + dt, dt)

# Resolvendo
Y = spi.odeint(EDOs, Y0, t)

# Separando as variaveis
x, th, vx, vth = Y.T

# Transformando as coordenadas
xm = L * np.sin(th) + x
ym = - L * np.cos(th)

# Plot
plt.figure()
plt.plot(t, x, 'r', label=r'$x(t)$')
plt.plot(t, th, 'b', label=r'$\theta(t)$')
plt.legend()
plt.show()


plt.figure()
plt.plot(t, xm, 'r', label=r'$xm(t)$')
plt.plot(t, ym, 'b', label=r'$ym(t)$')
plt.legend()
plt.show()
