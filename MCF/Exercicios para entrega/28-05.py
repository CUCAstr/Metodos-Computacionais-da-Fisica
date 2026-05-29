import numpy as np
import scipy.sparse as sps
import scipy.linalg as la
import matplotlib.pyplot as plt

# =========================================================
# Parâmetros Gerais e Discretização
# =========================================================
L = 1.0  # m
T = 1000.0  # N
N = 250  # discretização
x = np.linspace(0, L, N+1)
dx = L / N
x_int = x[1:-1] # nós internos (1 até N-1)

# Matriz M (diferenças finitas para derivada segunda)
# Diagonal principal = -2, diagonais vizinhas = 1
M = -sps.diags([np.ones(N-2), -2 * np.ones(N-1), np.ones(N-2)], [-1, 0, 1]).toarray()

# =========================================================
# EXERCÍCIO 1: Corda Homogênea
# =========================================================
rho_homog = 0.954  # kg/m
c2_homog = T / rho_homog  # (m/s)²

# Matriz de sobreposição constante para o caso homogêneo
S_homog = (dx**2 / c2_homog) * np.identity(N-1)

# Autovalores e autovetores
lbd_homog, X_homog = la.eigh(M, S_homog)
w_homog = np.sqrt(lbd_homog)
f_homog = w_homog / (2 * np.pi)

# =========================================================
# EXERCÍCIO 2: Corda Não Homogênea
# =========================================================
rho0 = 0.954  # kg/m
# δ = 0.5 g/m² = 0.5e-3 kg/m²
delta = 0.5e-3  # kg/m²

# Densidade variável dependente da posição x
# ρ = ρ0 + (x - L/2)*δ
rho_x = rho0 + (x_int - L/2) * delta

# Matriz de sobreposição diagonal usando densidade variável
# S_nn = dx² / c_n² = dx² * ρ(x_n) / T
S_nao_homog = np.diag((dx**2 * rho_x) / T)

# Autovalores e autovetores
lbd_nao_homog, X_nao_homog = la.eigh(M, S_nao_homog)
w_nao_homog = np.sqrt(lbd_nao_homog)
f_nao_homog = w_nao_homog / (2 * np.pi)

# =========================================================
# RESPOSTAS (a) e (b)
# =========================================================

# --- (a) Variação da Frequência Fundamental ---
var_f = f_nao_homog[0] - f_homog[0]

print("=== Resultados (a): Variação da Frequência Fundamental ===")
print(f"Frequência Ex1 (Homogênea):     {f_homog[0]:.7f} Hz")
print(f"Frequência Ex2 (Não Homogênea): {f_nao_homog[0]:.7f} Hz")
print(f"Variação (Ex2 - Ex1):           {var_f:.7e} Hz")

# --- (b) Comparação Gráfica do 1º Modo Normal ---
u_homog = np.zeros(N+1)
u_homog[1:-1] = X_homog[:, 0]  # 1º modo homogêneo

u_nao_homog = np.zeros(N+1)
u_nao_homog[1:-1] = X_nao_homog[:, 0] # 1º modo não homogêneo

# Normalização para comparação da forma de onda (amplitude máxima = 1)
u_homog = u_homog / np.max(np.abs(u_homog))
u_nao_homog = u_nao_homog / np.max(np.abs(u_nao_homog))

# Alinhamento de fase (garante que ambas as ondas comecem "para cima")
if u_homog[1] < 0:
    u_homog = -u_homog
if u_nao_homog[1] < 0:
    u_nao_homog = -u_nao_homog

plt.figure(figsize=(10, 6))
plt.plot(x, u_homog, label='Ex 1: Corda Homogênea', linestyle='--', color='blue')
plt.plot(x, u_nao_homog, label='Ex 2: Corda Não Homogênea', color='red', alpha=0.7)
plt.title('Comparação do Primeiro Modo Normal')
plt.xlabel('Posição x (m)')
plt.ylabel('Amplitude Normalizada')
plt.legend()
plt.grid(True)
plt.show()
