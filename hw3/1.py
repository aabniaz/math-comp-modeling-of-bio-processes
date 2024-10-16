# import numpy as np
# import matplotlib.pyplot as plt

# u = 1
# alpha = 1
# dx = 0.1
# dt = 0.00001
# eps = 1e-5
# L = 1.0
# Nx = int(L / dx) + 1
# Nt = 5000

# x = np.linspace(0, 1, Nx)
# T = np.zeros((Nt, Nx))

# T[0, :] = x*(x**2-3)
# T[:, 0] = 0
# T[:, -1] = 0

# CFL = (u * dt) / (dx ** 2)
# assert CFL <= 0.5, "CFL condition not satisfied"

# for n in range(Nt - 1):
#     for i in range(1, Nx - 1):
#         T[n + 1, i] = T[n, i] + dt * ((alpha**2 / dx**2) * (T[n, i + 1] - 2 * T[n, i] + T[n, i - 1]))

# x_values = np.linspace(0, 1, Nx)
# U_steady_state = T[-1, :]

# plt.plot(x_values, U_steady_state)
# plt.xlabel('x')
# plt.ylabel('U(t=steady state, x)')
# plt.show()








# import numpy as np
# import matplotlib.pyplot as plt

# alpha = 1
# dx = 0.1
# dt = 0.00001
# eps = 1e-5
# L = 1.0
# Nx = int(L / dx) + 1
# Nt = 5000

# x = np.linspace(0, 1, Nx)
# T = np.zeros((Nt, Nx))

# T[0, :] = x * (x**2 - 3)

# CFL = (alpha * dt) / (dx ** 2)
# assert CFL <= 0.5, "CFL condition not satisfied"

# for n in range(Nt - 1):
#     for i in range(1, Nx - 1):
#         T[n + 1, i] = T[n, i] + (alpha**2*dt/dx**2)* (T[n, i + 1] - 2 * T[n, i] + T[n, i - 1])

# # Neumann boundary conditions
# T[:, 0] = T[:, 1]
# T[:, -1] = T[:, -2]

# x_values = np.linspace(0, 1, Nx)
# U_steady_state = T[-1, :]

# plt.plot(x_values, U_steady_state)
# plt.xlabel('x')
# plt.ylabel('U(t=steady state, x)')
# plt.show()



# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters
# u = 1
# alpha = 1
# dx = 0.1
# dt = 0.00001
# eps = 1e-5
# L = 1.0
# Nx = int(L / dx) + 1
# Nt = 5000

# # Numerical solution
# x = np.linspace(0, 1, Nx)
# T = np.zeros((Nt, Nx))

# T[0, :] = x * (x**2 - 3)

# CFL = (u * dt) / (dx ** 2)
# assert CFL <= 0.5, "CFL condition not satisfied"

# for n in range(Nt - 1):
#     for i in range(1, Nx - 1):
#         T[n + 1, i] = T[n, i] + (dt*alpha**2 / dx**2) * (T[n, i + 1] - 2 * T[n, i] + T[n, i - 1])

# # Neumann boundary conditions
# T[:, 0] = T[:, 1]
# T[:, -1] = T[:, -2]

# # Analytical solution
# u_analytical = np.zeros(Nx)
# for n in range(1000):  
#     C = 2 * (-1)**n * (4 / ((2*n+1) * np.pi) + 48 / ((2*n+1)**3 * np.pi**3) + 96 / ((2*n+1)**4 * np.pi**4))
#     X = np.sin((2*n+1) * np.pi * np.linspace(0, 1, Nx) / 2)
#     T_analytical = np.exp(-((2*n+1)**2 * np.pi**2 * dt * Nt) / 4)
#     u_analytical += -C * X * T_analytical

# # Plotting
# plt.plot(x, T[-1, :], label='Numerical Solution (t=steady state)')
# plt.plot(x, u_analytical, label='Analytical Solution (t=steady state)', linestyle='--')

# plt.xlabel('x')
# plt.ylabel('Temperature (U)')
# plt.legend()
# plt.show()




import numpy as np
import matplotlib.pyplot as plt

def implicit(Nx, Nt, deltax, deltat):
    A = -1 / (deltax**2)
    B = 1 / deltat + 2 / (deltax**2)
    C = -1 / (deltax**2)

    alpha = np.zeros(Nx)
    beta = np.zeros(Nx)

    alpha[1] = 0
    beta[1] = 0

    # u^n_i (initial condition)
    u_n = np.linspace(0, 1, Nx)**3 - 3 * np.linspace(0, 1, Nx)

    U_implicit = np.zeros((Nt, Nx))
    U_implicit[0, :] = u_n
    
    for i in range(1, Nx-1):
        D = u_n[i] / deltat
        denominator = B + C * alpha[i]
        alpha[i+1] = -A / denominator
        beta[i+1] = (D - C * beta[i]) / denominator

    U_implicit[-1, :] = beta[Nx-1] / (1 - alpha[Nx-1])
    U = np.zeros(Nx)
    U[Nx-1] = U_implicit[-1, Nx-1]

    # U for interior points
    for i in range(Nx-1, 1, -1):
        U[i-1] = alpha[i] * U[i] + beta[i]

    return U

def analytical(x, t, N_terms):
    u_analytical = np.zeros_like(x)
    for n in range(1, N_terms + 1):
        C = 2 * (-1) ** (n - 1) * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
                    (2 * n - 1) ** 4 * np.pi ** 4))
        X = np.sin((2 * n - 1) * np.pi * x / 2)
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
        u_analytical += -C * X * T_analytical
    return u_analytical

def explicit(Nx, Nt, deltax, deltat, alpha):
    U_explicit = np.zeros((Nt, Nx))
    U_explicit[0, :] = np.linspace(0, 1, Nx)**3 - 3 * np.linspace(0, 1, Nx)  # Set initial condition

    for n in range(Nt - 1):
        for i in range(1, Nx - 1):
            U_explicit[n + 1, i] = U_explicit[n, i] + (alpha**2 * deltat / deltax**2) * (U_explicit[n, i+1] - 2 * U_explicit[n, i] + U_explicit[n, i-1])
        # Print intermediate values for debugging
            print(f'n={n}, i={i}, U_explicit={U_explicit[n + 1, i]}')


    return U_explicit

Nx = 100  
Nt = 500
L = 1.0
T = 0.5
alpha = 1
deltax = 0.1
deltat = 0.0001  
N_terms = 20

U_implicit = implicit(Nx, Nt, deltax, deltat)
U_explicit = explicit(Nx, Nt, deltax, deltat, alpha)
x_values = np.linspace(0, 1, Nx)
t_values = np.linspace(0, 0.5, Nt)
u_analytical = np.zeros((Nx, Nt))
for j, t in enumerate(t_values):
    u_analytical[:, j] = analytical(x_values, t, N_terms)

plt.figure(figsize=(12, 8))
plt.plot(x_values, U_implicit, label='Implicit', linestyle='--')
plt.plot(x_values, U_explicit[-1, :], label='Explicit', linestyle='--')
plt.plot(x_values, u_analytical[:, -1], label='Analytical', linestyle='-')
plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')
plt.legend()
plt.show()
print(U_explicit.shape)
