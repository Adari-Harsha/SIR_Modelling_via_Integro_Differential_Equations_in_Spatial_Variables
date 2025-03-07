import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Defining the Euler method
def euler_method(beta, gamma, tau, x_points, y_points, t_points, kernel_values):
    # Initialize arrays to store solutions
    z_solutions = np.zeros((len(t_points), len(x_points), len(y_points)))
    r_solutions = np.zeros((len(t_points), len(x_points), len(y_points)))

    # Setting our initial conditions
    center_x = len(x_points) // 2
    center_y = len(y_points) // 2
    #z_solutions[0, center_x, center_y] = 0.1
    z_solutions[0, center_x - 4:center_x + 6, center_y - 4:center_y + 6] = 0.1
    r_solutions[0, :, :] = 0

    # for loop
    for t_index in range(1, len(t_points)):
        # Update z using the discretized equation
        for i in range(len(x_points)):
            for j in range(len(y_points)):
                integrand = z_solutions[t_index - 1, :, :] * kernel_values[:, :, i, j]
                integral_term = np.trapz(np.trapz(integrand, x_points), y_points)

                z_solutions[t_index, i, j] = z_solutions[t_index - 1, i, j] + tau * (
                    beta * (1 - z_solutions[t_index - 1, i, j] - r_solutions[t_index - 1, i, j]) * integral_term - gamma * z_solutions[t_index - 1, i, j])

        # Updating r using the discretized equation
        r_solutions[t_index, :, :] = r_solutions[t_index - 1, :, :] + tau * gamma * z_solutions[t_index - 1, :, :]
        # Boundary Conditions
        z_solutions[t_index, [0, -1], :] = 0
        z_solutions[t_index, :, [0, -1]] = 0
        r_solutions[t_index, [0, -1], :] = 0
        r_solutions[t_index, :, [0, -1]] = 0


    return z_solutions, r_solutions

# Here we can tune our parameters
beta = 0.3
gamma = 0.1
tau = 0.5

# Defining our Spatial points (Grid size 30x30)
x_points = np.linspace(0, 1, 30)
y_points = np.linspace(0, 1, 30)

# Time points
t_points = np.arange(0, 200, tau)

# Initialize an array to store the kernel values
kernel_values = np.zeros((len(x_points), len(y_points), len(x_points), len(y_points)))

# Calculating our Kernel
for i in range(len(x_points)):
    for j in range(len(y_points)):
        for k in range(len(x_points)):
            for l in range(len(y_points)):
                kernel_values[i, j, k, l] = 500 * np.exp(-100 * np.sqrt((x_points[i] - x_points[k])**2 + (y_points[j] - y_points[l])**2))
                kernel_values[k, l, i, j] = kernel_values[i, j, k, l]


z_solutions, r_solutions = euler_method(beta, gamma, tau, x_points, y_points, t_points, kernel_values)

# Modify our graphs
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

def update_plot(frame):
    ax.clear()
    X, Y = np.meshgrid(x_points, y_points, indexing='ij')
    ax.plot_surface(X, Y, z_solutions[frame, :, :], cmap='viridis', rstride=1, cstride=1, alpha=0.8)
    ax.set_title(f'Time: {t_points[frame]:.1f}')
def init():
    ax.set_xlabel('Spatial Variable (x)')
    ax.set_ylabel('Spatial Variable (y)')
    ax.set_zlabel('Infections (z)')
    return ax

animation = FuncAnimation(fig, update_plot, frames=len(t_points), interval=100)
animation.save('infection_animation_2d.gif', writer='Pillow')
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(14, 6))

X, Y = np.meshgrid(x_points, y_points)
Z_final = z_solutions[-1, :, :]
R_final = r_solutions[-1, :, :]

# Plot for z
surf1 = ax1.plot_surface(X, Y, Z_final, cmap='viridis', edgecolor='none')
ax1.set_title('Final Infection Distribution')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('Infections (z)')
fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=5)

# Plot for r
surf2 = ax2.plot_surface(X, Y, R_final, cmap='plasma', edgecolor='none')
ax2.set_title('Final Recovery Distribution')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('Recoveries (r)')
fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=5)

plt.tight_layout()
plt.show()





