import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import matplotlib.pyplot as plt

T = 5000 # N
rho0 = 4 # kg/m
g = 9.81 # m/s^2
L = 8 # m

m1, a1, x1 = 40, 0.3, 3.5 # kg, m, m
m2, a2, x2 = 55, 0.27, 4.1 # kg, m, m

'''
rho1 = (m1 / 2 * a1) * (np.exp(-((np.abs(x-x1) / a1)**2))) # kg/m
rho2 = (m2 / 2 * a2) * (np.exp(-((np.abs(x-x2) / a2)**2))) # kg/m
'''

N = 250 # 250 nos
x = np.linspace(0, L, N+1)


trunc1 = np.where(np.abs(x-x1) <= a1)
trunc2 = np.where(np.abs(x-x2) <= a2)   #trunc retorna true para os indices onde a condição é satisfeita e false para os outros indices. Assim, g1 e g2 serão preenchidos apenas nos indices onde a condição é satisfeita, ou seja, onde x está dentro do intervalo definido por a1 e a2 em torno de x1 e x2, respectivamente.

g1 = np.zeros(N+1)
g2 = np.zeros(N+1)
g1[trunc1] = (m1 / 2 * a1) * (np.exp(-((np.abs(x[trunc1]-x1) / a1)**2))) # kg/m
g2[trunc2] = (m2 / 2 * a2) * (np.exp(-((np.abs(x[trunc2]-x2) / a2)**2))) # kg/m


Dx = L / N # Delta x
M = sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray() # Matriz Diagonal

C = g * Dx**2 / (T / (rho0 + g1 + g2)[1:-1]) # g/c^2

U = np.zeros(N+1) # deslocamentos de todos os nos da corda
U[1:-1] = la.solve(M, C) # deslocamentos dos nos internos



plt.figure()
plt.plot(x, U)
plt.show()