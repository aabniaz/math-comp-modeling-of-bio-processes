import numpy as np
import matplotlib.pyplot as plt

dx, dy, dt = 0.027, 0.027, 0.0001
Re = 1
rho = 1
nu = 1 / Re
eps, eps_P = 1e-6, 1e-6
stop_it, stop_it_P = 1e4, 1e6
alpha = 1

x = np.arange(0, 1 + dx, dx)
y = np.arange(0, 1 + dy, dy)
N = len(x)
M = len(y)

u_old = np.zeros((M, N))
u_star = np.zeros((M, N))
u_new = np.zeros((M, N))

v_old = np.zeros((M, N))
v_star = np.zeros((M, N))
v_new = np.zeros((M, N))

p_old = np.zeros((M, N))
p_new = np.zeros((M, N))

T = np.zeros((M, N))

u_old[:, :] = 0
v_old[:, :] = 0
p_old[:, :] = 0

m1 = int(0.3/dy)
m2 = int(0.6/dy)

# Inlet for u and v
u_old[-1, m2:] = 0
v_old[-1, m2:] = 1

# Outlet for u and v
u_old[0, m2:] = 0
u_old[m1:m2, 0] = 0
v_old[0, m2:] = 0
v_old[m1:m2, 0] = 0

# Walls for u and v
u_old[:, -1] = 0
u_old[:m1, 0] = 0
u_old[m2:, 0] = 0
u_old[0, :m2] = 0
u_old[-1, :m2] = 0

v_old[:, -1] = 0
v_old[:m1, 0] = 0
v_old[m2:, 0] = 0
v_old[0, :m2] = 0
v_old[-1, :m2] = 0

T_out = -10
T_hot = 50

T[0, :m1] = T_hot # heat source
T[-1, m2:] = T_out # inlet
T[0, m2:] = T[1, m2:] # outlet1
T[m1:m2, 0] = T[m1:m2, 1] # outlet 2
# walls
T[:, -1] = (T[:, -2] + T_out)/ 2 
T[:m1, 0] = (T[:m1, 1] + T_out) / 2
T[m2:, 0] = (T[m2:, 1] + T_out) / 2
T[0, m1:m2] = (T[1, m1:m2] + T_out) / 2
T[-1, :m2] = (T[-2, :m2] + T_out) / 2

maximum = 1
iteration = 0
while maximum > eps and iteration < stop_it:
        #u_star, v_star
        u_star[1:-1, 1:-1] = u_old[1:-1, 1:-1] - dt*(u_old[1:-1, 1:-1]*(u_old[1:-1, 2:] - u_old[1:-1, 0:-2])/(2*dx) + v_old[1:-1, 1:-1]*(u_old[2:, 1:-1] - u_old[0:-2, 1:-1])/(2*dy)- nu*( (u_old[1:-1, 2:] - 2*u_old[1:-1, 1:-1] + u_old[1:-1, 0:-2])/dx**2 + (u_old[2:, 1:-1] - 2*u_old[1:-1, 1:-1] + u_old[0:-2, 1:-1])/dy**2)) 

        v_star[1:-1, 1:-1] = v_old[1:-1, 1:-1] - dt*(u_old[1:-1, 1:-1]*(v_old[1:-1, 2:] - v_old[1:-1, 0:-2])/(2*dx) + v_old[1:-1, 1:-1]*(v_old[2:, 1:-1] - v_old[0:-2, 1:-1])/(2*dy) - nu*((v_old[1:-1, 2:] - 2*v_old[1:-1, 1:-1] + v_old[1:-1, 0:-2])/dx**2 + (v_old[2:, 1:-1] - 2*v_old[1:-1, 1:-1] + v_old[0:-2, 1:-1])/dy**2)) 

        u_star[-1, m2:] = 0
        u_star[0, m2:] = 0
        u_star[m1:m2, 0] = 0
        u_star[:, -1] = 0
        u_star[:m1, 0] = 0
        u_star[m2:, 0] = 0
        u_star[0, :m2] = 0
        u_star[-1, :m2] = 0

        v_star[-1, m2:] = 1
        v_star[0, m2:] = 0
        v_star[m1:m2, 0] = 0
        v_star[:, -1] = 0
        v_star[:m1, 0] = 0
        v_star[m2:, 0] = 0
        v_star[0, :m2] = 0
        v_star[-1, :m2] = 0

        #pressure poisson eq
        maximum_P = 1
        iteration_P = 0
        while maximum_P > eps_P and iteration_P < stop_it_P:
                p_new[1:-1, 1:-1] = (dy**2*(p_old[1:-1, 2:] + p_old[1:-1, 0:-2]) + dx**2*(p_old[2:, 1:-1] + p_old[0:-2, 1:-1]))/(2*(dx**2 + dy**2)) - dx**2*dy**2*rho/(2*dt*(dx**2 + dy**2)) * ((u_star[1:-1, 2:] - u_star[1:-1, 0:-2])/(2*dx) + (v_star[2:, 1:-1] - v_star[0:-2, 1:-1])/(2*dy)) 

                #boundary cond for p
                p_new[:, -1] = p_new[:, -2]
                p_new[:m1, 0] = p_new[:m1, 1]
                p_new[m2:, 0] = p_new[m2:, 1]
                p_new[0, :m2] = p_new[1, :m2]
                p_new[-1, :m2] = p_new[-2, :m2]
                # inlet
                p_new[-1, m2:] = 1
                # outlet
                p_new[0, m2:] = 0
                p_new[m1:m2, 0] = 0

                maximum_P = np.max(np.abs(p_new - p_old))
                p_old[:, :] = p_new[:, :]
                iteration_P += 1

        # corrector
        u_new[1:-1, 1:-1] = u_star[1:-1, 1:-1] - dt/rho * (p_new[1:-1, 2:] - p_new[1:-1, 0:-2])/(2*dx)   
        v_new[1:-1, 1:-1] = v_star[1:-1, 1:-1] - dt/rho * (p_new[2:, 1:-1] - p_new[0:-2, 1:-1])/(2*dy) 

        u_new[-1, m2:] = 0
        u_new[0, m2:] = 0
        u_new[m1:m2, 0] = 0
        u_new[:, -1] = 0
        u_new[:m1, 0] = 0
        u_new[m2:, 0] = 0
        u_new[0, :m2] = 0
        u_new[-1, :m2] = 0

        v_new[-1, m2:] = 1
        v_new[0, m2:] = 0
        v_new[m1:m2, 0] = 0
        v_new[:, -1] = 0
        v_new[:m1, 0] = 0
        v_new[m2:, 0] = 0
        v_new[0, :m2] = 0
        v_new[-1, :m2] = 0

        u_old[:, :] = u_new[:, :]
        v_old[:, :] = v_new[:, :]

        T_new = np.zeros_like(T)
        T_new[1:-1, 1:-1] = dt * (alpha**2*((T[1:-1, 2:] - 2*T[1:-1, 1:-1] + T[1:-1, :-2])/dx**2 + (T[2:, 1:-1] - 2*T[1:-1, 1:-1] + T[:-2, 1:-1])/dy**2) - u_new[1:-1, 1:-1]*((T[1:-1, 2:] - T[1:-1, 1:-1])/dx) - v_new[1:-1, 1:-1]*((T[2:, 1:-1] - T[1:-1, 1:-1])/dy)) + T[1:-1, 1:-1]

        T_new[0, :m1] = T_hot # heat source
        T_new[-1, m2:] = T_out # inlet
        T_new[0, m2:] = T_new[1, m2:] # outlet1
        T_new[m1:m2, 0] = T_new[m1:m2, 1] # outlet 2
        # walls
        T_new[:, -1] = (T_new[:, -2] + T_out)/ 2 
        T_new[:m1, 0] = (T_new[:m1, 1] + T_out) / 2
        T_new[m2:, 0] = (T_new[m2:, 1] + T_out) / 2
        T_new[0, m1:m2] = (T_new[1, m1:m2] + T_out) / 2
        T_new[-1, :m2] = (T_new[-2, :m2] + T_out) / 2

        max_t = np.max(np.abs(T - T_new))
        T[:, :] = T_new[:, :]

        iteration += 1

X, Y = np.meshgrid(x, y)

# print(iter)

# fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# axs[0, 0].contourf(X, Y, u_new, origin='lower', extent=(0, 1, 0, 1))
# axs[0, 0].set_title('u')

# axs[0, 1].contourf(X, Y, v_new, origin='lower', extent=(0, 1, 0, 1))
# axs[0, 1].set_title('v')

# axs[1, 0].contourf(X, Y, p_new, origin='lower', extent=(0, 1, 0, 1))
# axs[1, 0].set_title('p')

# speed = np.sqrt(u_new**2 + v_new**2)
# axs[1, 1].contourf(X, Y, speed, origin='lower', extent=(0, 1, 0, 1))
# axs[1, 1].set_title('u + v')

# plt.tight_layout()
# plt.show()

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, T, cmap='viridis')
# plt.colorbar(label='T')
plt.title('T')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()