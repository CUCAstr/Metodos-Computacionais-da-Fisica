import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def resolver_catenaria(xA, yA, xB, yB, L, rho, g, guess=[2.0, -1.0, 70.0]):
    """
    Resolve os parâmetros da curva da catenária exata.
    """
    def equations(vars):
        x0, y0, T0 = vars
        a = T0 / (rho * g)
        
        # Condições de contorno dos pontos A e B
        eq1 = y0 + a * np.cosh((xA - x0) / a) - yA
        eq2 = y0 + a * np.cosh((xB - x0) / a) - yB
        
        # Condição do comprimento total
        eq3 = a * (np.sinh((xB - x0) / a) - np.sinh((xA - x0) / a)) - L
        return [eq1, eq2, eq3]
        
    return fsolve(equations, guess)

def resolver_parabola(xA, yA, xB, yB, L, rho, g, guess=[1.8, -1.1, 64.0]):
    """
    Resolve os parâmetros da aproximação parabólica usando integral analítica fechada.
    """
    def equations(vars):
        C1, C2, T0 = vars
        a = (rho * g) / T0
        
        # Condições de contorno dos pontos A e B
        eq1 = C1 + C2 * xA + (a / 2.0) * xA**2 - yA
        eq2 = C1 + C2 * xB + (a / 2.0) * xB**2 - yB
        
        # Integral do comprimento da parábola calculada analiticamente
        def F(x):
            u = C2 + a * x
            return (u * np.sqrt(1.0 + u**2) + np.arcsinh(u)) / (2.0 * a)
            
        eq3 = (F(xB) - F(xA)) - L
        return [eq1, eq2, eq3]
        
    return fsolve(equations, guess)

# Execução com os parâmetros de exemplo da Figura 4
xA, yA = 1.0, 1.0
xB, yB = 4.0, 2.0
L = 3.5
rho = 4.0
g = 9.81

# Resolução dos sistemas
x0_c, y0_c, T0_c = resolver_catenaria(xA, yA, xB, yB, L, rho, g)
C1_p, C2_p, T0_p = resolver_parabola(xA, yA, xB, yB, L, rho, g)

# Saída dos dados
print("=== Catenária (Solução Exata) ===")
print(f"x0 = {x0_c:.5f} m")
print(f"y0 = {y0_c:.5f} m")
print(f"T0 = {T0_c:.5f} N\n")

print("=== Parábola (Aproximação) ===")
print(f"C1 = {C1_p:.5f} m")
print(f"C2 = {C2_p:.5f} m")
print(f"T0 = {T0_p:.5f} N")

# --- Plotagem dos Resultados ---
x_plot = np.linspace(xA, xB, 100)
a_c = T0_c / (rho * g)
y_cat = y0_c + a_c * np.cosh((x_plot - x0_c) / a_c)

a_p = (rho * g) / T0_p
y_par = C1_p + C2_p * x_plot + (a_p / 2.0) * x_plot**2

plt.figure(figsize=(8, 5))
plt.plot(x_plot, y_cat, 'b-', label='Catenária (Exata)', linewidth=2)
plt.plot(x_plot, y_par, 'r--', label='Parábola (Aproximação)', linewidth=1.5)
plt.scatter([xA, xB], [yA, yB], color='black', zorder=5, label='Pontos A e B')

plt.title('Comparação: Catenária vs Aproximação Parabólica')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.show()