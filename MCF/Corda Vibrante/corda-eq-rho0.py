import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import matplotlib.pyplot as plt

T = 5000  # N
rho0 = 4  # kg/m
g = 9.81  # m/s²
L = 8     # m

c2 = T / rho0  # c² (m²/s²)
gc2 = g / c2  # g / c²  (1/m)
N = 250  # 250 nós
Dx = L / N  # Delta x
M = sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray()

C = gc2 * np.ones(N-1) * Dx**2

U = np.zeros(N+1)  # deslocamentos de todos os nós da corda
U[1:-1] = la.solve(M, C)  # deslocamentos dos nós internos

x = np.linspace(0, L, N+1)

plt.figure()
plt.plot(x, U)
plt.show()
