import numpy as np
import matplotlib.pyplot as plt

dx = 0.01
dy = 0.01
dt = 0.00001
Re = 1
rho = 1
nu = 1/Re
alpha = 0.8

e = 1e-6
stop_iteration = 10000
stop_iteration_P = 10000

x = np.arange(0, 1 + dx, dx)
y = np.arange(0, 1 + dy, dy)
N = len(x)
M = len(y)
print(N, M)

T = np.zeros((M, N))

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

def bound_t(t):
    m1 = int(0.4/dy)
    m2 = int(0.6/dy)
    
    T_out = -10
    T_hot = 50
    # T out != T hot
    
    # heat source
    t[-1, m1:m2] = T_hot
    
    # walls
    t[-1, :m1] = (t[-2, :m1] + T_out)/2
    t[-1, m2:] = (t[-2, m2:] + T_out)/2
    t[0, :m1] = (t[1, :m1] + T_out)/2
    t[0, m2:] = (t[1, m2:] + T_out)/2
    t[:m1, 0] = (t[:m1, 1] + T_out)/2
    t[m2:, 0] = (t[m2:, 1] + T_out)/2
    t[:m1, -1] = (t[:m1, -2] + T_out)/2
    t[m2:, -1] = (t[m2, -2] + T_out)/2
    
    # inlet
    t[m1:m2, -1] = T_out
    
    # outlet
    t[m1:m2, 0] = t[m1:m2, 1]
    t[0, m1:m2] = t[1, m1:m2]


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
bound_t(T)

max = 1
max_t = 1
iter = 0
while max > e  and max_t > e and iter < stop_iteration:
    
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
    
    # temperature     
    T_new = np.zeros_like(T)
    T_new[1:-1, 1:-1] = dt * (alpha**2*((T[1:-1, 2:] - 2*T[1:-1, 1:-1] + T[1:-1, :-2])/dx**2 + (T[2:, 1:-1] - 2*T[1:-1, 1:-1] + T[:-2, 1:-1])/dy**2) - U_new[1:-1, 1:-1]*((T[1:-1, 2:] - T[1:-1, 1:-1])/dx) - V_new[1:-1, 1:-1]*((T[2:, 1:-1] - T[1:-1, 1:-1])/dy)) + T[1:-1, 1:-1]
    bound_t(T_new)
    
    max_t = np.max(np.abs(T - T_new))
    T[:, :] = T_new[:, :]
        
    print(max_t)
    
    iter += 1

X, Y = np.meshgrid(x, y)
print(iter)

# надо добавить графики для T, U, V, P