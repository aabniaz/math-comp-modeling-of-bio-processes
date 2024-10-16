import numpy as np
import matplotlib.pyplot as plt

dx = 0.01
dy = 0.01
dt = 0.0001
Re = 20
rho = 1
nu = 1/Re

e = 1e-6
stop_iteration = 10000
stop_iteration_P = 1e6

x = np.arange(0, 1 + dx, dx)
y = np.arange(0, 1 + dy, dy)
N = len(x)
M = len(y)
print(N, M)

U_old = np.zeros((M, N))
U_s = np.zeros((M, N))
U_new = np.zeros((M, N))

V_old = np.zeros((M, N))
V_s = np.zeros((M, N))
V_new = np.zeros((M, N))

P_old = np.zeros((M, N))
P_new = np.zeros((M, N))

# at t=0:
U_old[:, :] = 0
V_old[:, :] = 0
P_old[:, :] = 0

def bound_u(u):
    m1 = int(0.4/dy)
    m2 = int(0.6/dy)

    # walls
    u[-1, :] = 0
    u[0, :m1] = 0
    u[0, m2:] = 0
    u[:m1, 0] = 0
    u[m2:, 0] = 0
    u[:m1, -1] = 0
    u[m2:, -1] = 0
    
    # inlet
    u[m1:m2, -1] = -1
    
    # outlet
    u[m1:m2, 0] = u[m1:m2, 1]
    u[0, m1:m2] = u[1, m1:m2]

def bound_v(v):
    m1 = int(0.4/dy)
    m2 = int(0.6/dy)
    
    # walls
    v[-1, :] = 0
    v[0, :m1] = 0
    v[0, m2:] = 0
    v[:m1, 0] = 0
    v[m2:, 0] = 0
    v[:m1, -1] = 0
    v[m2:, -1] = 0
    
    # inlet
    v[m1:m2, -1] = 0
    
    # outlet
    v[m1:m2, 0] = v[m1:m2, 1]
    v[0, m1:m2] = v[1, m1:m2]

def bound_p(p):
    m1 = int(0.4/dy)
    m2 = int(0.6/dy)
    
    # walls
    p[-1, :] = p[-2, :]
    p[0, :m1] = p[1, :m1]
    p[0, m2:] = p[1, m2:]
    p[:m1, 0] = p[:m1, 1]
    p[m2:, 0] = p[m2:, 1]
    p[:m1, -1] = p[:m1, -2]
    p[m2:, -1] = p[m2:, -2]
    
    # inlet
    p[m1:m2, -1] = -1
    
    # outlet
    p[m1:m2, 0] = 0
    p[0, m1:m2] = 0
    
bound_u(U_old)
bound_v(V_old)
bound_p(P_old)

max = 1
iter = 0
while max > e and iter < stop_iteration:
    
    #step 1: Burgers eq
    U_s[1:-1, 1:-1] = U_old[1:-1, 1:-1] - dt*(
        U_old[1:-1, 1:-1]*(U_old[1:-1, 2:] - U_old[1:-1, 0:-2])/(2*dx) \
        + V_old[1:-1, 1:-1]*(U_old[2:, 1:-1] - U_old[0:-2, 1:-1])/(2*dy) \
        - nu*(
            (U_old[1:-1, 2:] - 2*U_old[1:-1, 1:-1] + U_old[1:-1, 0:-2])/dx**2 \
            + (U_old[2:, 1:-1] - 2*U_old[1:-1, 1:-1] + U_old[0:-2, 1:-1])/dy**2))
    
    V_s[1:-1, 1:-1] = V_old[1:-1, 1:-1] - dt*(
        U_old[1:-1, 1:-1]*(V_old[1:-1, 2:] - V_old[1:-1, 0:-2])/(2*dx) \
        + V_old[1:-1, 1:-1]*(V_old[2:, 1:-1] - V_old[0:-2, 1:-1])/(2*dy) \
        - nu*(
            (V_old[1:-1, 2:] - 2*V_old[1:-1, 1:-1] + V_old[1:-1, 0:-2])/dx**2 \
            + (V_old[2:, 1:-1] - 2*V_old[1:-1, 1:-1] + V_old[0:-2, 1:-1])/dy**2))
    
    bound_u(U_s)
    bound_v(V_s)

    
    max_P = 1
    iter_P = 0
    
    # step 2: Poissons eq
    while max_P > e and iter_P < stop_iteration_P:

        P_new[1:-1, 1:-1] = (dy**2*(P_old[1:-1, 2:] + P_old[1:-1, 0:-2]) \
            + dx**2*(P_old[2:, 1:-1] + P_old[0:-2, 1:-1]))/(2*(dx**2 + dy**2)) \
                - dx**2*dy**2*rho/(2*dt*(dx**2 + dy**2)) \
                * ((U_s[1:-1, 2:] - U_s[1:-1, 0:-2])/(2*dx) \
                + (V_s[2:, 1:-1] - V_s[0:-2, 1:-1])/(2*dy))

        bound_p(P_new)
        
        max_P = np.max(np.abs(P_new - P_old))
#         print(max_P, iter_P)
        
        P_old[:, :] = P_new[:, :]
        
        iter_P += 1
    
    # step 3: Corrector
    U_new[1:-1, 1:-1] = U_s[1:-1, 1:-1] - dt/rho * (P_new[1:-1, 2:] - P_new[1:-1, 0:-2])/(2*dx)  
    V_new[1:-1, 1:-1] = V_s[1:-1, 1:-1] - dt/rho * (P_new[2:, 1:-1] - P_new[0:-2, 1:-1])/(2*dy)

    bound_u(U_new)
    bound_v(V_new)
        
    max = np.max(np.abs(U_new - U_old))
    U_old[:, :] = U_new[:, :]
    V_old[:, :] = V_new[:, :]
    # print(max)
    
    iter += 1

X, Y = np.meshgrid(x, y)

print(iter)

import matplotlib.pyplot as plt

# Create figure and axes
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Plot u
axs[0, 0].contourf(X, Y, U_new, origin='lower', extent=(0, 1, 0, 1))
axs[0, 0].set_title('U')

# Plot v
axs[0, 1].contourf(X, Y, V_new, origin='lower', extent=(0, 1, 0, 1))
axs[0, 1].set_title('V')

# Plot p
axs[1, 0].contourf(X, Y, P_new, origin='lower', extent=(0, 1, 0, 1))
axs[1, 0].set_title('P')

# Plot |u| + |v|
speed = np.sqrt(U_new**2 + V_new**2)
axs[1, 1].contourf(X, Y, speed, origin='lower', extent=(0, 1, 0, 1))
axs[1, 1].set_title('U + V')

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()