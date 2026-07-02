import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

def minus_V(vars):
    R, H = vars
    return - np.pi * R**2 * H

V_max = opt.minimize(fun = minus_V,
                     x0 = [5,5],
                     bounds = [(5,20),(0,20)],
                     constraints = {'type': 'ineq', 'fun': lambda x: 900 - (2*np.pi * x[0] * x[1])})

print(V_max)



# ==========================================
# CÓDIGO DO GRÁFICO DE CURVAS DE NÍVEL
# ==========================================

# 1. Criando a "Malha" (Grid)
# Vamos gerar valores de Raio (R) e Altura (H) de 0 até 25 para dar um bom panorama
vetor_R = np.linspace(0, 25, 200)
vetor_H = np.linspace(0, 25, 200)
malha_R, malha_H = np.meshgrid(vetor_R, vetor_H)

# 2. Calculando Volume e Área para cada ponto da malha
malha_Volume = np.pi * (malha_R**2) * malha_H
malha_Area = 2 * np.pi * malha_R * malha_H

plt.figure(figsize=(8, 6))

# 3. Desenhando as Curvas de Nível do Volume (montanhas)
# levels=20 significa desenhar 20 linhas de "altitude" diferentes
contorno_vol = plt.contour(malha_R, malha_H, malha_Volume, levels=20, cmap='viridis')
plt.colorbar(contorno_vol, label='Volume (cm³)') # Adiciona a barra de cores ao lado

# 4. Desenhando a restrição da Área e suas próprias curvas de nível
# Desenha 10 curvas de área, tracejadas, mais finas e levemente transparentes
plt.contour(malha_R, malha_H, malha_Area, 10, colors='red', linestyles='dashed', alpha=0.4)

# Desenha a linha EXATA da nossa restrição (Área = 900) bem forte e contínua
plt.contour(malha_R, malha_H, malha_Area, levels=[900], colors='red', linewidths=3)

# 5. Desenhando os Limites da "Caixa" (Bounds)
plt.axvline(x=5, color='black', linestyle='--', label='R mín (5 cm)')
plt.axvline(x=20, color='black', linestyle='--', label='R máx (20 cm)')
plt.axhline(y=20, color='black', linestyle='--', label='H máx (20 cm)')

# 6. Marcando a sua resposta (O Ponto Ótimo)
R_otimo = V_max.x[0]
H_otimo = V_max.x[1]
plt.plot(R_otimo, H_otimo, 'r*', markersize=15, label='Volume Máximo')

# 7. Acabamentos estéticos
plt.title('Curvas de Nível: Maximização do Volume da Caneca')
plt.xlabel('Raio R (cm)')
plt.ylabel('Altura H (cm)')
plt.xlim(0, 25)
plt.ylim(0, 25)
plt.legend(loc='upper left')

# Mostra o gráfico na tela
plt.show()