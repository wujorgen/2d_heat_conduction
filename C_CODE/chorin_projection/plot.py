import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt("test_pressure_sol.csv", delimiter=',')
u = np.loadtxt("test_u_sol.csv", delimiter=',')
v = np.loadtxt("test_v_sol.csv", delimiter=',')


N_POINTS = 41
DOMAIN_SIZE = 1.0
x = np.linspace(0.0, DOMAIN_SIZE, p.shape[1])
y = np.linspace(0.0, DOMAIN_SIZE, p.shape[0])
X, Y = np.meshgrid(x, y)

plt.contourf(X, Y, p)
plt.colorbar()
plt.streamplot(X, Y, u, v, color="black")
plt.title("Lid Driven Cavity")
plt.show()

