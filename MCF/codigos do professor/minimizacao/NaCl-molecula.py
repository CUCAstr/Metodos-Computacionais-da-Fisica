import numpy as np
import scipy.optimize as spo

# Dados
V0 = 1.09e3  # eV
r0 = .330  # Angstrom
e2 = 14.4  # eV Ang

def V(r):
    return -e2 / r + V0 * np.exp(-r / r0)

sol = spo.minimize(V, 3)
print(sol)
