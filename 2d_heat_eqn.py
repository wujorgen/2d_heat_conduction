import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

plate_length = 100
plate_width = 50
max_iter_time = 1000

alpha = 2.0
delta_x = 1
delta_y = 1

delta_t = (delta_x ** 2) / (4 * alpha)
gamma = (alpha * delta_t) / (delta_x ** 2)

u = np.empty( (max_iter_time, plate_length, plate_width) )

u_initial = 0.0

u_top = 100.0
u_left = 0.0
u_bottom = 0.0
u_right = 0.0

u.fill(u_initial)

u[:, (plate_length-1):, :] = u_top
u[:, :, :1] = u_left
u[:, :1, 1:] = u_bottom
u[:, :, (plate_length-1):] = u_right

def calculate(u):
    for k in range(0, max_iter_time - 1, 1):
        print(f"On step {k} of {max_iter_time}")
        for i in range(1, plate_length - 1, delta_x):
            for j in range(1, plate_width - 1, delta_y):
                u[k + 1, i, j] = alpha * delta_t  \
                                * ( \
                                    (u[k, i + 1, j] + u[k, i - 1, j] - 2*u[k, i, j]) / (delta_x**2) + \
                                    (u[k, i, j + 1] + u[k, i, j - 1] - 2*u[k, i, j]) / (delta_y**2) \
                                ) \
                                + u[k,i,j]
                #u[k + 1, i, j] = gamma * (u[k][i+1][j] + u[k][i-1][j] + u[k][i][j+1] + u[k][i][j-1] - 4*u[k][i][j]) + u[k][i][j]

    return u

def plotheatmap(u_k, k):
    # Clear the current plot figure
    plt.clf()
    plt.title(f"Temperature at t = {k*delta_t:.3f} unit time")
    plt.xlabel("x")
    plt.ylabel("y")
    
    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=0, vmax=100)
    plt.colorbar()
    
    return plt

def animate(k):
    plotheatmap(u[k], k)


u = calculate(u)
anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=max_iter_time, repeat=False)
anim.save("heat_equation_solution.gif")