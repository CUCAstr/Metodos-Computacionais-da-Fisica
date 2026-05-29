import numpy as np
import scipy.sparse as sps
import scipy.linalg as la
import matplotlib.pyplot as plt

# Parâmetros Físicos
L = 8.0  # m
T = 5000.0  # N
rho0 = 4.0  # kg/m

# Dados dos amigos (massa, centro, raio)
m1, x1, a1 = 40.0, 3.5, 0.3
m2, x2, a2 = 55.0, 4.1, 0.27

# Malha de discretização
N = 250
x = np.linspace(0, L, N+1)
dx = L / N

# Cálculo da densidade ponto a ponto (apenas nós internos)
x_int = x[1:-1]
rho_n = np.ones(N-1) * rho0

# Adição da densidade da Gaussiana truncada para Amigo 1
idx1 = np.abs(x_int - x1) <= a1
rho_n[idx1] += (m1 / (2 * a1)) * np.exp(-((x_int[idx1] - x1)**2) / a1**2)

# Adição da densidade da Gaussiana truncada para Amigo 2
idx2 = np.abs(x_int - x2) <= a2
rho_n[idx2] += (m2 / (2 * a2)) * np.exp(-((x_int[idx2] - x2)**2) / a2**2)

# Construção da Matriz M exatamente como em corda_homog.py
M = -sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray()

# Construção da Matriz S (Matriz de sobreposição)
# Sabendo que S da aula usa (dx**2) / c2 e c2 = T / rho(x)
# Logo: S = (dx**2) * rho(x) / T
S = np.diag(dx**2 * rho_n / T)

# Problema Generalizado de Autovalores e Autovetores
lbd, X = la.eigh(M, S)

# Cálculo das Frequências
w = np.sqrt(lbd)
freqs = w / (2 * np.pi) # Conversão de rad/s para Hz

print("5 primeiras frequências de vibração (Hz):")
print(freqs[:5])

# Plotagem dos autoestados (5 primeiros modos)
plt.figure(figsize=(10, 6))
for i in range(5):
    u = np.zeros(N+1)
    u[1:-1] = X[:, i]  # i-ésimo modo
    
    # Normalização gráfica do vetor para melhor visualização (opcional)
    u = u / np.max(np.abs(u))
    
    plt.plot(x, u, label=f"Modo {i+1} ({freqs[i]:.2f} Hz)")

plt.title("Os 5 Primeiros Modos Normais da Corda")
plt.xlabel("Posição x (m)")
plt.ylabel("Amplitude Normalizada")
plt.legend()
plt.grid(True)
plt.show()