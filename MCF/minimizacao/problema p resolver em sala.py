import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

#variaveis
a, b = 0.15, 0.2
c = 1-(a+b)
p = 10
C = 1
K = 5

# Funcao
def lucro(vars):
    L, M = vars
    return -(p * C * K**a * L**b * M**c - L - M)


sol = opt.minimize(lucro, [1,1])
print(sol)