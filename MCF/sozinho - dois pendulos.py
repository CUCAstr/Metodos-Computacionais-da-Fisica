# %%
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# %% [markdown]
# Coordenadas das massas:
# 
# $x_1 = l \sin(\theta_1)$
#  
# $y_1 = -l \cos(\theta_1)$
#  
# $x_2 = l_0 + l \sin(\theta_2)$
#   
# $y_2 = -l \cos(\theta_2)$
# 
# Derivando as coordenadas:
#
# $\dot{x}_1 = l \cos(\theta_1) \dot{\theta}_1$
#
# $\dot{y}_1 = l \sin(\theta_1) \dot{\theta}_1$
#
# $\dot{x}_2 = l \cos(\theta_2) \dot{\theta}_2$
#
# $\dot{y}_2 = l \sin(\theta_2) \dot{\theta}_2$
# %%

l = 1    # Comprimento dos pêndulos (m)
g = 10   # Aceleração da gravidade (m/s^2)
m1 = 1   # Massa do pêndulo 1 (kg)
m2 = 1   # Massa do pêndulo 2 (kg)
k = 20   # Constante da mola (N/m)
l0 = 2   # Comprimento de repouso da mola (m)

theta1_0 = 0.5  # Pêndulo 1 inclinado (rad)
theta2_0 = 0.0  # Pêndulo 2 em repouso
Vtheta1_0 = 0.0  # Velocidade angular inicial do pêndulo 1 (rad/s)
Vtheta2_0 = 0.0  # Velocidade angular inicial do pêndulo 2 (rad/s)

Y = [theta1_0, Vtheta1_0, theta2_0, Vtheta2_0]

x1 = l * np.sin(theta1_0)
y1 = -l * np.cos(theta1_0)
x2 = l0 + l * np.sin(theta2_0)
y2 = -l * np.cos(theta2_0)

def pendulos_acoplados(Y, t):
    theta1, theta2, Vtheta1, Vtheta2 = Y
    Dtheta1 = Vtheta1
    Dtheta2 = Vtheta2
    DVtheta1 = -(m * g * l * np.sin(theta1) - k * np.sqrt((x2 - x1)**2 + (y2 - y1)**2) - l0 )/ (m * l**2)


