# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters
# L = 1.0  # Length of the rod
# T = 0.5  # Total simulation time
# Nx = 100  # Number of spatial grid points
# Nt = 500  # Number of time steps

# # Spatial and temporal discretization
# dx = L / (Nx - 1)
# dt = T / Nt

# # Initial condition
# def initial_condition(x):
#     return x**3 - 3*x

# # Initialize solution matrix
# u = np.zeros((Nx, Nt+1))

# # Set initial condition
# u[:, 0] = initial_condition(np.linspace(0, L, Nx))

# # Finite difference coefficients
# A = -1 / dx**2
# B = 1 / dt + 2 / dx**2
# C = -1 / dx**2

# # Time-stepping loop using implicit method and Thomas algorithm
# for n in range(0, Nt):
#     # Prepare the tridiagonal system
#     D = u[:, n] / dt
#     D[0] = 0  # Boundary condition
#     D[-1] = 0  # Boundary condition

#     # Thomas algorithm
#     alpha = np.zeros(Nx-2)
#     beta = np.zeros(Nx-2)

#     for i in range(1, Nx-1):
#         denominator = B + C * alpha[i-1]
#         alpha[i-1] = -A / denominator
#         beta[i-1] = (D[i] - C * beta[i-1]) / denominator

#     # Backward substitution
#     for i in range(Nx-3, -1, -1):
#         u[i+1, n+1] = alpha[i] * u[i+2, n+1] + beta[i]

# # Plot the results
# x_values = np.linspace(0, L, Nx)
# t_values = np.linspace(0, T, Nt+1)

# X, T = np.meshgrid(x_values, t_values)

# fig = plt.figure(figsize=(10, 6))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, T, u.T, cmap='viridis')
# ax.set_xlabel('Position (x)')
# ax.set_ylabel('Time (t)')
# ax.set_zlabel('Temperature (u)')
# ax.set_title('Numerical Solution of 1D Heat Equation with Initial Condition x^3 - 3x (Implicit Method)')
# plt.show()












import numpy as np
import matplotlib.pyplot as plt

L = 1.0  
T = 0.5  
Nx = 100 
Nt = 500  
N_terms = 20  

dx = L / (Nx - 1)
dt = T / Nt

def initial_condition(x):
    return x**3 - 3*x

# def analytical_solution(x, t, N_terms):
#     u_analytical = np.zeros_like(x)
#     for n in range(N_terms):
#         term = (8 / ((2*n+1) * np.pi)) * np.cos((2*n+1) * np.pi / 2) * (1 + 12 / ((2*n+1)**2 * np.pi**2))
#         term *= np.sin((2*n+1) * np.pi * x / 2) * np.exp(-((2*n+1)**2 * np.pi**2 / 4) * t)
#         u_analytical += term
#     return u_analytical


def analytical_solution(x, t, N_terms):
    u_analytical = np.zeros_like(x)
    for n in range(N_terms):
        C = 2 * (-1)**n * (4 / ((2*n+1) * np.pi) + 48 / ((2*n+1)**3 * np.pi**3) + 96 / ((2*n+1)**4 * np.pi**4))
        X = np.sin((2*n+1) * np.pi * x / 2)
        T = np.exp(-((2*n+1)**2 * np.pi**2 * t) / 4)
        u_analytical += C * X * T
    return u_analytical
u = np.zeros((Nx, Nt+1))

u[:, 0] = initial_condition(np.linspace(0, L, Nx))

A = -1 / dx**2
B = 1 / dt + 2 / dx**2
C = -1 / dx**2

# for n in range(0, Nt):
#     D = u[:, n] / dt
#     D[0] = 0  
#     D[-1] = 0 

#     alpha = np.zeros(Nx-2)
#     beta = np.zeros(Nx-2)

#     for i in range(1, Nx-1):
#         denominator = B + C * alpha[i-1]
#         alpha[i-1] = -A / denominator
#         beta[i-1] = (D[i] - C * beta[i-1]) / denominator

#     for i in range(Nx-3, -1, -1):
#         u[i+1, n+1] = alpha[i] * u[i+2, n+1] + beta[i]


for n in range(0, Nt):
    # Prepare the tridiagonal system
    D = u[1:Nx-1, n] / dt  # Modified to exclude boundary points
    D[0] -= A * u[0, n+1]  # Adjust for the boundary condition
    D[-1] -= C * u[Nx-1, n+1]  # Adjust for the boundary condition

    alpha = np.zeros(Nx-2)
    beta = np.zeros(Nx-2)

    # Thomas algorithm
    for i in range(1, Nx-2):
        denominator = B + C * alpha[i-1]
        alpha[i] = -A / denominator
        beta[i] = (D[i] - C * beta[i-1]) / denominator

    # Backward substitution
    for i in range(Nx-3, -1, -1):
        u[i+1, n+1] = alpha[i] * u[i+2, n+1] + beta[i]


x_values = np.linspace(0, L, Nx)
t_values = np.linspace(0, T, Nt+1)

X, T = np.meshgrid(x_values, t_values)
fig = plt.figure(figsize=(12, 8))

ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, T, u.T, cmap='viridis')
ax1.set_title('Numerical Solution')
ax1.set_xlabel('Position (x)')
ax1.set_ylabel('Time (t)')
ax1.set_zlabel('Temperature (u)')

u_analytical = np.zeros((Nx, Nt+1))
for j, t in enumerate(t_values):
    u_analytical[:, j] = analytical_solution(x_values, t, N_terms)

ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, T, u_analytical.T, cmap='viridis')
ax2.set_title('Analytical Solution')
ax2.set_xlabel('Position (x)')
ax2.set_ylabel('Time (t)')
ax2.set_zlabel('Temperature (u)')

plt.show()











# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters
# L = 1.0  # Length of the rod
# T = 0.5  # Total simulation time
# Nx = 100  # Number of spatial grid points
# Nt = 500  # Number of time steps
# N_terms = 20  # Number of terms in the Fourier series

# # Spatial and temporal discretization
# dx = L / (Nx - 1)
# dt = T / Nt

# # Initial condition
# def initial_condition(x):
#     return x**3 - 3*x

# # Analytical solution at a specific time t
# def analytical_solution(x, t, N_terms):
#     u_analytical = np.zeros_like(x)
#     for n in range(N_terms):
#         term = (8 / ((2*n+1) * np.pi)) * np.cos((2*n+1) * np.pi / 2) * (1 + 12 / ((2*n+1)**2 * np.pi**2))
#         term *= np.sin((2*n+1) * np.pi * x / 2) * np.exp(-((2*n+1)**2 * np.pi**2 / 4) * t)
#         u_analytical += term

#     return u_analytical



# # Initialize solution matrix
# u = np.zeros((Nx, Nt+1))

# # Set initial condition
# u[:, 0] = initial_condition(np.linspace(0, L, Nx))

# # Finite difference coefficients
# A = -1 / dx**2
# B = 1 / dt + 2 / dx**2
# C = -1 / dx**2

# # Time-stepping loop using implicit method and Thomas algorithm
# for n in range(0, Nt):
#     # Prepare the tridiagonal system
#     D = u[:, n] / dt
#     D[0] = 0  # Boundary condition
#     D[-1] = 0  # Boundary condition

#     # Thomas algorithm
#     alpha = np.zeros(Nx-2)
#     beta = np.zeros(Nx-2)

#     for i in range(1, Nx-1):
#         denominator = B + C * alpha[i-1]
#         alpha[i-1] = -A / denominator
#         beta[i-1] = (D[i] - C * beta[i-1]) / denominator

#     # Backward substitution
#     for i in range(Nx-3, -1, -1):
#         u[i+1, n+1] = alpha[i] * u[i+2, n+1] + beta[i]

# # Plot the results
# x_values = np.linspace(0, L, Nx)
# t_values = np.linspace(0, T, Nt+1)

# X, T = np.meshgrid(x_values, t_values)

# fig = plt.figure(figsize=(12, 8))

# # Plot Numerical Solution
# ax1 = fig.add_subplot(121, projection='3d')
# ax1.plot_surface(X, T, u.T, cmap='viridis')
# ax1.set_title('Numerical Solution')
# ax1.set_xlabel('Position (x)')
# ax1.set_ylabel('Time (t)')
# ax1.set_zlabel('Temperature (u)')

# # Plot Analytical Solution at a specific time (e.g., T/2 for comparison)
# t_analytical = T / 2
# u_analytical = analytical_solution(x_values, t_analytical, N_terms)

# ax2 = fig.add_subplot(122)
# ax2.plot(x_values, u_analytical, label='Analytical Solution')
# ax2.plot(x_values, u[:, -1], label='Numerical Solution', linestyle='--')
# ax2.set_title('Comparison at t = {:.3f}'.format(t_analytical))
# ax2.set_xlabel('Position (x)')
# ax2.set_ylabel('Temperature (u)')
# ax2.legend()

# plt.show()
