import numpy as np
import scipy.optimize as spo
import scipy.integrate as spi
import matplotlib.pyplot as plt

#  Dados
xa, ya = 1, 1  # Ponto A
xb, yb = 5, 2  # Ponto B

# Solução exata (círculo)

def circulo(x, y, x0, r):
    return (x - x0) ** 2 + y ** 2 - r ** 2

def ajusta_circulo(pars):
    eq1 = circulo(xa, ya, *pars)
    eq2 = circulo(xb, yb, *pars)
    return [eq1, eq2]

sol = spo.root(ajusta_circulo, [1, 1])
x0, r = sol.x

N = 100
x = np.linspace(xa, xb, N+1)
y = np.array([spo.newton(lambda y: circulo(xn, y, x0, r), x0) for xn in x])

TempoExato = r * spi.simpson(1 / y**2, x)

# solução numérica

def TempoAprox(yt):
    yt[0] = ya
    yt[-1] = yb
    dx = np.diff(x)
    dy = np.diff(yt)
    ds = np.sqrt(dx ** 2 + dy ** 2)
    vy = yt[:-1]
    dt = ds / vy
    return dt.sum()

# y0 = np.interp(xc, [xa, xb], [ya, yb])
sol = spo.minimize(TempoAprox, y)
yap = sol.x
taprox = sol.fun

# Gráficos

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.plot(x, y, '-', label=f'Tempo exato: {TempoExato:.3f}')
ax.plot(x, yap, '.', label=f'Tempo aproximado: {taprox:.3f}')
ax.legend()
plt.tight_layout()
plt.show()
