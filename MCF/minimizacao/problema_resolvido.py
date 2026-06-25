import numpy as np
import scipy.optimize as opt

# Variáveis (parâmetros do problema Cobb-Douglas)
a, b = 0.15, 0.2
c = 1 - (a + b)
p = 10
C = 1
K = 5

# Função de Lucro
# Minimizamos o negativo do lucro para encontrar o máximo
def lucro(vars):
    L, M = vars
    
    # Proteção extra contra valores não-positivos (evita números complexos)
    if L <= 0 or M <= 0:
        return 1e18
    
    # Receita (p * Produção) - Custos (L + M)
    f_producao = C * (K**a) * (L**b) * (M**c)
    receita = p * f_producao
    custo = L + M
    
    return -(receita - custo)

# 1. Definimos limites (bounds) para garantir que L e M sejam sempre positivos
# (1e-6 é um valor pequeno próximo de zero)
limites = [(1e-6, None), (1e-6, None)]

# 2. Chute inicial (ajustado para a escala real do problema, aproximadamente 10^5)
chute_inicial = [100000, 100000]

# 3. Execução da otimização usando o método L-BFGS-B
# Adicionamos 'tol' para ajudar o algoritmo a convergir na escala correta
sol = opt.minimize(lucro, chute_inicial, method='L-BFGS-B', bounds=limites, tol=1e-8)

# Exibição dos resultados
print("--- Resultado da Otimização ---")
if sol.success:
    print(f"Sucesso: {sol.success}")
    print(f"L otimizado: {sol.x[0]:.2f}")
    print(f"M otimizado: {sol.x[1]:.2f}")
    print(f"Lucro máximo: {-sol.fun:.2f}")
else:
    print("A otimização falhou.")
    print(f"Mensagem: {sol.message}")

print("\n--- Objeto de Solução Completo ---")
print(sol)
