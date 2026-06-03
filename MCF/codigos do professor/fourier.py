'''
    Módulo fourier
    Cálculos de séries e transformadas de Fourier
    Consulte a ajuda de cada função para mais informações
'''

import numpy as np
import scipy.integrate as spi

class serie:
    '''
        Classe fourier.serie
        Cálculo numérico de uma série de Fourier

        $    f(t) = a_0 + \\sum_{n=1}^nmax a_n \\cos (2 \\pi n t / T)
                        + \\sum_{n=1}^nmax b_n \\sin n (2 \\pi n t / T)
        $

        Uso:

        Criação de Objeto:
            serie(f, t1, t2, nmax)
                f: a função python a ser expandida (deve depender de
                   um único argumento)
                t1: limite de integração inferior (float, padrão: -pi)
                t2: limite de integração superior (float, padrão: +pi)
                nmax: máximo n nas somas em seno e cosseno (int, padrão: 5)
        Métodos:
            .coeficientes(f, t1, t2, nmax)
                retorna os coeficientes da série de Fourier como
                    self.a e self.b (floats)
            .calc(t)
                calcula a série para um certo valor de t (float ou numpy array)
        Propriedades:
            Da criação do objeto:
                .f: a função fornecida pelo usuário na criação do objeto (function)
                .T: semiperíodo da série
                .nmax: máximo n nas somas em seno e cosseno (int)
                .a: coeficientes dos cossenos (array) - inclui a0
                .b: coeficientes dos senos (array)
    '''

    def __init__(self, f, t1=-np.pi, t2=np.pi, nmax=5):
        self.f = f
        self.T = .5 * (t2 - t1)
        self.nmax = nmax
        self.coeficientes(f, t1, t2, nmax)
        
    def coeficientes(self, f, t1, t2, nmax):
        self.a = np.zeros(nmax+1)
        self.b = np.zeros(nmax+1)
        self.a[0] = self.a0_coef(f, t1, t2)
        for n in range(1, nmax+1):
            self.a[n] = self.an_coef(f, n, t1, t2)
            self.b[n] = self.bn_coef(f, n, t1, t2)

    def a0_coef(self, f, t1, t2):
        v, err = spi.quad(f, t1, t2)
        T = .5 * (t2 - t1)
        return v / T

    def an_coef(self, f, n, t1, t2):

        def fc(t):
            return f(t) * np.cos(np.pi * n * t / T)

        T = .5 * (t2 - t1)
        v, err = spi.quad(fc, t1, t2)
        return v / T

    def bn_coef(self, f, n, t1, t2):

        def fs(t):
            return f(t) * np.sin(np.pi * n * t / T)

        T = .5 * (t2 - t1)
        v, err = spi.quad(fs, t1, t2)
        return v / T

    def calc(self, t):
        k = np.pi * t / self.T
        # conterá o valor da série
        if type(t) == float or type(t) == int:
            serie = .5 * self.a[0]
        else:  # supondo t um array ou lista
            serie = np.zeros_like(t)
            serie[:] = .5 * self.a[0]
        for n in range(1, self.nmax+1):
            ang = k * n
            serie += self.a[n] * np.cos(ang)
            serie += self.b[n] * np.sin(ang)
        return serie
