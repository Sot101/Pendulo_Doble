import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation

# Longitudes de las varillas del péndulo (m) y masas de los péndulos (kg)
L1, L2 = 1, 1
m1, m2 = 1, 1
# Aceleración gravitacional (m/s^2)
g = 9.81

def deriv(t, y, L1, L2, m1, m2):
    """Devuelve las primeras derivadas de y = theta1, z1, theta2, z2."""
    theta1, z1, theta2, z2 = y

    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1dot = z1
    z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    theta2dot = z2
    z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
             m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    return [theta1dot, z1dot, theta2dot, z2dot]

def calc_E(y):
    """Devuelve la energía total del sistema."""
    th1, th1d, th2, th2d = y.T
    V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
    T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
            2*L1*L2*th1d*th2d*np.cos(th1-th2))
    return T + V

# Tiempo máximo, intervalos de tiempo y la cuadrícula de tiempo (todo en segundos)
tmax, dt = 30, 0.001  # Reduce el tamaño del paso
t = np.arange(0, tmax+dt, dt)
# Condiciones iniciales: theta1, dtheta1/dt, theta2, dtheta2/dt
y0 = np.array([3*np.pi/7, 0, 3*np.pi/4, 0])

# Realizar la integración numérica de las ecuaciones de movimiento
sol = solve_ivp(deriv, [0, tmax], y0, args=(L1, L2, m1, m2), t_eval=t, method='RK45', rtol=1e-12, atol=1e-12)
y = sol.y.T

# Verificar que el cálculo conserva la energía total dentro de cierta tolerancia
EDRIFT = 0.1  # Permitir un mayor desvío de energía
# Energía total a partir de las condiciones iniciales
E = calc_E(np.array([y0]))
if np.max(np.sum(np.abs(calc_E(y) - E))) > EDRIFT:
    raise ValueError('Se excedió el máximo desvío de energía de {}.'.format(EDRIFT))

# Desempaquetar z y theta como una función del tiempo
theta1, theta2 = y[:,0], y[:,2]

# Convertir a coordenadas cartesianas las posiciones de los dos péndulos
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# Radio del círculo que representa los péndulos
r = 0.05
# Traza de la posición del péndulo m2 durante los últimos trail_secs segundos
trail_secs = 1
# Esto corresponde a max_trail puntos de tiempo
max_trail = int(trail_secs / dt)

fig, ax = plt.subplots(figsize=(8, 6))

def init():
    ax.set_xlim(-L1-L2-r, L1+L2+r)
    ax.set_ylim(-L1-L2-r, L1+L2+r)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    return []

def update(i):
    ax.clear()
    ax.set_xlim(-L1-L2-r, L1+L2+r)
    ax.set_ylim(-L1-L2-r, L1+L2+r)
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    
    # Las varillas del péndulo
    ax.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], lw=2, c='k')
    # Círculos que representan el punto de anclaje de la varilla 1 y los péndulos 1 y 2
    c0 = Circle((0, 0), r/2, fc='k', zorder=10)
    c1 = Circle((x1[i], y1[i]), r, fc='b', ec='b', zorder=10)
    c2 = Circle((x2[i], y2[i]), r, fc='r', ec='r', zorder=10)
    ax.add_patch(c0)
    ax.add_patch(c1)
    ax.add_patch(c2)
    
    # La traza se dividirá en ns segmentos y se graficará como una línea que se desvanece
    ns = 20
    s = max_trail // ns

    for j in range(ns):
        imin = i - (ns-j)*s
        if imin < 0:
            continue
        imax = imin + s + 1
        # El desvanecimiento se ve mejor si se cuadran las longitudes fraccionarias a lo largo de la traza
        alpha = (j/ns)**2
        ax.plot(x2[imin:imax], y2[imin:imax], c='r', solid_capstyle='butt',
                lw=2, alpha=alpha)
    return []

fps = 30  # Ajustar la tasa de fotogramas
ani = FuncAnimation(fig, update, frames=range(0, len(t), int(1/fps/dt)), init_func=init, blit=True)

plt.show()
