import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt

# Dados
V0 = 1.09e3  # eV
r0 = .330  # Angstrom
e2 = 14.4  # eV Ang
n, m = 5, 5

# Criando os íons

class ion:
    # Classe ion: cria o objeto
    # ion, de carga .carga na posição .pos
    def __init__(self, carga):
        self.carga = carga
        self.pos = [0, 0, 0]  # ao criar o ion, coloca-o na origem

Ions = []  # lista com os íons
for i in range(n):
    Ions.append(ion(1))
for j in range(m):
    Ions.append(ion(-1))

# "chute" inicial: íons em posições aleatórias
pos0 = np.random.rand((n + m) * 3)  # será posteriormente atribuída aos íons

# Potencial de dois íons
def V(I1, I2):
    s = I1.carga * I2.carga  # -1: atração; 1: repulsão
    V0a = V0 if I1.carga - I2.carga != 0 else 0
    r = np.linalg.norm(I1.pos - I2.pos)
    return s * e2 / r + V0a * np.exp(-r / r0)

# Potencial total do cluster
def Vcluster(pos):

    position = pos.reshape(n + m, 3)
    for i in range(len(Ions)):
        Ions[i].pos = position[i]

    Vc = 0
    for i in range(len(Ions)):
        for j in range(i+1, len(Ions)):
            Vc += V(Ions[i], Ions[j])

    return Vc

# Encontrando as posições de equilíbrio
sol = spo.minimize(Vcluster, pos0)
print(sol)

# reformatando as posições para um formato mais amigável
x, y, z = sol.x.reshape(n + m, 3).T

# Gráficos
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
colors = [i.carga for i in Ions]
ax.scatter(x, y, z, c=colors, cmap='jet')

plt.show()
