import numpy.random as rd
import numpy as np
import math

n_tentativas = 10**5
n_pessoas = 16
soma = 0

for i in range(n_tentativas):
    
    aniv = []
    
    # Geração de aniversários
    for j in range(n_pessoas):
        aniv.append(rd.randint(1,366))
        j += 1
    
    # Check de pelo menos 1 data repetida
    if len(aniv) > len(np.unique(aniv)):
        soma += 1
    
    i += 1

# Probabilidade teórica
P_teorico = 1 - (math.factorial(365)/((math.factorial(365-n_pessoas))*(365**n_pessoas)))

# Probabilidade calculada
P_calculado = soma/n_tentativas

# Calculo do erro
erro = 1 - P_calculado/P_teorico


print(f"Pessoas: {n_pessoas}\nIterações: {n_tentativas}")
print(f"Teórico: {P_teorico}")
print(f"Calculado: {P_calculado}")
print(f"Erro: {erro}%")