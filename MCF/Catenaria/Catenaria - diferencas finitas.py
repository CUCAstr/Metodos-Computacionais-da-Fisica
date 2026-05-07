import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import matplotlib.pyplot as plt

T = 5000 # N
rho0 = 4 # kg/m
g = 9.81 # m/s^2
L = 8 # m

c2 = T / rho0 # c^2 (m/s^2)
gc2 = g / c2 # g/c^2

N= 250 # 250 nos
Dx = L / N # Delta x
M = sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray() # Matriz Diagonal

C = gc2 * np.ones(N-1) * Dx**2

U = np.zeros(N+1) # deslocamentos de todos os nos da corda
U[1:-1] = la.solve(M, C) # deslocamentos dos nos internos

plt.figure()
plt.plot(U)
plt.show()



