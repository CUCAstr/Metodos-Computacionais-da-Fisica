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

# solução numérica - método de Ritz

deg = 4  # grau do polinômio na solução teste

def sol_teste(x, a):
    p = a[0]
    for n in range(1, deg+1):
        p += a[n] * x ** (n)
    return p

def df(f, x, dx=1e-4):  # derivada numérica
    return (f(x+dx) - f(x)) / dx

def Obj_f(a):
    y = sol_teste(x, a)
    dy = df(lambda x: sol_teste(x, a), x)
    f = np.sqrt(1 + dy ** 2) / y
    return spi.simpson(f, x)

eq_consa = {'type': 'eq', 'fun': lambda a: sol_teste(xa, a) - ya}
eq_consb = {'type': 'eq', 'fun': lambda a: sol_teste(xb, a) - yb}
constraints = [eq_consa, eq_consb]

sol = spo.minimize(Obj_f, x0=np.ones(deg+1), constraints=constraints)
yap = sol_teste(x, sol.x)
taprox = sol.fun
print(sol.x)

# Gráficos

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.plot(x, y, '-', label=f'Tempo exato: {TempoExato:.3f}')
ax.plot(x, yap, '.', label=f'Tempo aproximado: {taprox:.3f}')
ax.legend()
plt.tight_layout()
plt.show()

