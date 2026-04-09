import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as pla
import scipy.integrate as spi

L0 = 3    #m
L  = 4    #m
g  = 9.81 #m/s
k  = 30   #N/m  


def pendulos_ligados(Y,t):
    tht1, tht2, vt1, vt2 = Y
    
    df_dtht1 = -2 * L * (L0 + L * np.sin(tht2) - L + np.sin(tht1)) * np.cos(tht1) + 2 * L * (L*np.cos(tht2)-L*np.cos(tht1))*np.sin(tht1)
    df_dtht2 =  2 * L * (L0 + L * np.sin(tht2) - L + np.sin(tht1)) * np.cos(tht2) - 2 * L * (L*np.cos(tht2)-L*np.cos(tht1))*np.sin(tht2)
    eps = np.sqrt((L0 + L * np.sin(tht2) - L * np.sin(tht1))*2 + (L * np.cos(tht2) - L * np.cos(tht1))*2) - L0
    deps_dtht1 = df_dtht1/(2*(eps+L0))
    deps_dtht2 = df_dtht2/(2*(eps+L0))
  
    dtht1 = vt1
    dtht2 = vt2
    dvt1  = -(g*np.sin(tht1)-k*eps*deps_dtht1)
    dvt2  = -(g*np.sin(tht2)-k*eps*deps_dtht2)
    return [dtht1, dtht2, dvt1, dvt2]

tht10, tht20, vt10, vt20 = np.pi/12, np.pi/9, 0, 0
Y0 = [tht10, tht20, vt10, vt20]

tmin, tmax, h = 0, 20, 0.01
t = np.arange(tmin, tmax+h, h)

Y = spi.odeint(pendulos_ligados, Y0, t)

tht1, tht2, vt1, vt2 = Y.T
x1 = L*np.sin(tht1)
y1 = - L*np.cos(tht1)

x2 = L0 + L * np.sin(tht2)
y2 = - L*np.cos(tht2)

plt.figure()
plt.plot(x1, y1, label='pendulo 1')
plt.plot(x2, y2, label='pendulo 2')
plt.legend()
plt.show()