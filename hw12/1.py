import numpy as np
import matplotlib.pyplot as plt

dx, dy, dt = 0.027, 0.027, 0.0001
Re = 1
rho = 1
nu = 1 / Re
eps, eps_P = 1e-6, 1e-6
speed_U, speed_P = 1, 1
stop_it, stop_it_P = 1e4,1e6

x = np.arange(0, 1 + dx, dx)
y = np.arange(0, 1 + dy, dy)
N = len(x)
M = len(y)
K = 100  
inlet_end = int(0.2 / dy)
outlet_start = int(0.8 / dy)

u = np.zeros((M, N))
u_star = np.zeros((M, N))
u_new = np.zeros((M, N))

v = np.zeros((M, N))
v_star = np.zeros((M, N))
v_new = np.zeros((M, N))

p = np.zeros((M, N))
p_new = np.zeros((M, N))

u[outlet_start:, 0] = speed_U
v[outlet_start:, 0] = 0
p[outlet_start:, 0] = speed_P

# walls
u[0:outlet_start, 0] = 0
u[inlet_end:outlet_start, N - 1] = 0
u[0, :] = 0
u[M - 1, :] = 0

v[0:outlet_start, 0] = 0
v[inlet_end:outlet_start, N - 1] = 0
v[0, :] = 0
v[M - 1, :] = 0

p[0:outlet_start, 0] = p[0:outlet_start, 1]
p[inlet_end:outlet_start, N - 1] = p[inlet_end:outlet_start, N - 2]
p[0, :] = p[1, :]
p[M - 1, :] = p[M - 2, :]

# outlet
u[0:inlet_end, N - 1] = u[0:inlet_end, N - 2]
u[outlet_start:, N - 1] = u[outlet_start:, N - 2]
v[0:inlet_end, N - 1] = 0
v[outlet_start:, N - 1] = 0
p[0:inlet_end, N - 1] = 0
p[outlet_start:, N - 1] = 0

maximum = 1
iteration = 0
while maximum > eps and iteration < stop_it:
    #u_star, v_star
    u_star[1:-1, 1:-1] = u[1:-1, 1:-1] - dt * (
            u[1:-1, 1:-1] * (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)
            + v[1:-1, 1:-1] * (u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy)
            - nu * (
            (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2]) / dx ** 2
            + (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1]) / dy ** 2
        )
    )

    v_star[1:-1, 1:-1] = v[1:-1, 1:-1] - dt * (
            u[1:-1, 1:-1] * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)
            + v[1:-1, 1:-1] * (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)
            - nu * (
            (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]) / dx ** 2
            + (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1]) / dy ** 2
        )
    )

    #pressure poisson eq
    maximum_P = 1
    iteration_P = 0
    while maximum_P > eps_P and iteration_P < stop_it_P:
        p_new[1:-1, 1:-1] = (dy ** 2 * (p[1:-1, 2:] + p[1:-1, 0:-2])
                            + dx ** 2 * (p[2:, 1:-1] + p[0:-2, 1:-1])) / (
                                   2 * (dx ** 2 + dy ** 2)) \
                           - dx ** 2 * dy ** 2 * rho / (
                                   2 * dt * (dx ** 2 + dy ** 2)) \
                           * ((u_star[1:-1, 2:] - u_star[1:-1, 0:-2]) / (2 * dx)
                              + (v_star[2:, 1:-1] - v_star[0:-2, 1:-1]) / (2 * dy))

        #boundary cond for p
        p_new[0:outlet_start, 0] = p_new[0:outlet_start, 1]
        p_new[outlet_start:, 0] = speed_P
        p_new[0:inlet_end, N - 1] = 0
        p_new[inlet_end:outlet_start, N - 1] = p_new[inlet_end:outlet_start, N - 2]
        p_new[outlet_start:, N - 1] = 0
        p_new[0, :] = p_new[1, :]
        p_new[M - 1, :] = p_new[M - 2, :]
        maximum_P = np.max(np.abs(p_new - p))
        p[:, :] = p_new[:, :]
        iteration_P += 1

    u_new[1:-1, 1:-1] = u_star[1:-1, 1:-1] - dt / rho * (
            p_new[1:-1, 2:] - p_new[1:-1, 0:-2]) / (2 * dx)
    v_new[1:-1, 1:-1] = v_star[1:-1, 1:-1] - dt / rho * (
            p_new[2:, 1:-1] - p_new[0:-2, 1:-1]) / (2 * dy)

    u_new[0:outlet_start, 0] = 0
    u_new[outlet_start:, 0] = 1
    u_new[0:inlet_end, N - 1] = u_new[0:inlet_end, N - 2]
    u_new[inlet_end:outlet_start, N - 1] = 0
    u_new[outlet_start:, N - 1] = u_new[outlet_start:, N - 2]
    u_new[0, :] = 0
    u_new[M - 1, :] = 0

    v_new[0:outlet_start, 0] = 0
    v_new[outlet_start:, 0] = 0
    v_new[0:inlet_end, N - 1] = v_new[0:inlet_end, N - 2]
    v_new[inlet_end:outlet_start, N - 1] = 0
    v_new[outlet_start:, N - 1] = v_new[outlet_start:, N - 2]
    v_new[0, :] = 0
    v_new[M - 1, :] = 0

    u[:, :] = u_new[:, :]
    v[:, :] = v_new[:, :]

    iteration += 1

plt.figure(figsize=(12, 4))

plt.subplot(131)
plt.contourf(x, y, u_new, cmap='viridis', levels=20)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('u')

plt.subplot(132)
plt.contourf(x, y, v_new, cmap='viridis', levels=20)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('v')

# plt.subplot(133)
# plt.contourf(x, y, p_new, cmap='viridis', levels=20)
# plt.colorbar()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('p')

plt.tight_layout()
plt.show()