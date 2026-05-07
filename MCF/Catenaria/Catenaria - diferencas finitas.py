import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import matplotlib.pyplot as plt

T = 5000 # N
rho0 = 4 # kg/m
g = 9.81 # m/s^2
L = 8 # m

m1, a1, x1 = 40, 0.3, 3.5
m2, a2, x2 = 55, 0.27, 7.1

N = 250 # 250 nos
x = np.linspace(0, L, N+1)

trunc1 = np.where(np.abs(x-x1) <= a1)
trunc2 = np.where(np.abs(x-x2) <= a2)

g1 = np.zeros(N+1)
g2 = np.zeros(N+1)
g1[trunc1] = .5 * m1 / a1 * np.exp(-((x[trunc1] - x1) / a1)**2)
g2[trunc2] = .5 * m2 / a2 * np.exp(-((x[trunc2] - x2) / a2)**2)

M = sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray() # Matriz Diagonal

dens = rho0 + g1 + g2
c2 = T / dens[1:-1] # c^2 (m/s^2)
Dx = L / N # Delta x
C = g * Dx**2 / c2

U = np.zeros(N+1) # deslocamentos de todos os nos da corda
U[1:-1] = la.solve(M, C) # deslocamentos dos nos internos


plt.figure()
plt.plot(x, U)
plt.show()



