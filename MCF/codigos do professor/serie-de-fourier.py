import numpy as np
import matplotlib.pyplot as plt
import fourier

'''
Cálculo da série de Fourier de uma função f(t) no intervalo [0, 8]
para um certo nmax
usando integração numérica
'''

# Para ver a ajuda do módulo e das classes, descomente as linhas abaixo:
# help(fourier)
# help(fourier.serie)
# (ou abra o arquivo fourier.py)

# definindo a função a ser expandida em série de Fourier:
def f(t):  # a função deve depender de um único argumento
    return t ** 2

# GRÁFICO

plt.figure()
plt.xlabel(f'$t$')
plt.ylabel(f'$f(t)$')
plt.axvline(0, color='k')
plt.axhline(0, color='k')

t = np.arange(-20, 20, .01)
y = f(t)
plt.plot(t, y, label=f'$f(t)$')


t1, t2 = 0, 8
nmax = 10
sf = fourier.serie(f, t1=t1, t2=t2, nmax=nmax)
serie = sf.calc(t)

plt.plot(t, serie, label=f'$s(t)$, $n_{{max}}={nmax}$')

plt.axvline(t1, color='k', linestyle='--')
plt.axvline(t2, color='k', linestyle='--')

plt.legend(loc='best')
plt.show()
