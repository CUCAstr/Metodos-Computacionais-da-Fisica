import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt

# dados
a = 1.3
w0 = 1.5
g = 9.81
l = 0.7
m = 0.5
th0 = np.pi/6

def f_acel(t, th):
    # equacao do pendulo no vagem
    return (a*w0**2*np.cos(w0*t)*np.cos(th) - g*np.sin(th))/l

# item a - rk4
tmax = 10.0
Nt = 1000
t = np.linspace(0, tmax, Nt)
dt = t[1]-t[0]

theta = np.zeros(Nt)
omega = np.zeros(Nt)
theta[0] = th0
omega[0] = 0

for i in range(Nt-1):
    # passo rk4
    k1th = omega[i]
    k1om = f_acel(t[i], theta[i])
    
    k2th = omega[i] + 0.5*dt*k1om
    k2om = f_acel(t[i] + 0.5*dt, theta[i] + 0.5*dt*k1th)
    
    k3th = omega[i] + 0.5*dt*k2om
    k3om = f_acel(t[i] + 0.5*dt, theta[i] + 0.5*dt*k2th)
    
    k4th = omega[i] + dt*k3om
    k4om = f_acel(t[i] + dt, theta[i] + dt*k3th)
    
    theta[i+1] = theta[i] + (dt/6)*(k1th + 2*k2th + 2*k3th + k4th)
    omega[i+1] = omega[i] + (dt/6)*(k1om + 2*k2om + 2*k3om + k4om)

plt.figure()
plt.plot(t, theta)
plt.title('theta(t)')
plt.show()

# item b - fft
T2 = 200.0
N2 = 2**13
t2 = np.linspace(0, T2, N2, endpoint=False)
dt2 = t2[1]-t2[0]

th2 = np.zeros(N2)
om2 = np.zeros(N2)
th2[0] = th0
om2[0] = 0

for i in range(N2-1):
    k1th = om2[i]
    k1om = f_acel(t2[i], th2[i])
    k2th = om2[i] + 0.5*dt2*k1om
    k2om = f_acel(t2[i]+0.5*dt2, th2[i]+0.5*dt2*k1th)
    k3th = om2[i] + 0.5*dt2*k2om
    k3om = f_acel(t2[i]+0.5*dt2, th2[i]+0.5*dt2*k2th)
    k4th = om2[i] + dt2*k3om
    k4om = f_acel(t2[i]+dt2, th2[i]+dt2*k3th)
    th2[i+1] = th2[i] + (dt2/6)*(k1th + 2*k2th + 2*k3th + k4th)
    om2[i+1] = om2[i] + (dt2/6)*(k1om + 2*k2om + 2*k3om + k4om)

trans = fft.fft(th2)
w_v = 2*np.pi*fft.fftfreq(N2, dt2)
espectro = np.abs(trans)**2

wp = np.sqrt(g/l)

plt.figure()
mask = (w_v > 0) & (w_v < 10)
plt.plot(w_v[mask], espectro[mask])
plt.axvline(w0, color='r', ls='--')
plt.axvline(wp, color='b', ls='--')
plt.show()

print("w0:", w0)
print("wp:", wp)
