# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters
# L = 1.0  # Length of the rod
# T = 0.5  # Total simulation time
# Nx = 100  # Number of spatial grid points
# N_terms = 20  # Number of terms in the Fourier series

# # Spatial discretization
# x_values = np.linspace(0, L, Nx)

# # Analytical solution at a specific time t
# def analytical_solution(x, t, N_terms):
#     u_analytical = np.zeros_like(x)
#     for n in range(N_terms):
#         term = (8 / ((2*n+1) * np.pi)) * np.cos((2*n+1) * np.pi / 2) * (1 + 12 / ((2*n+1)**2 * np.pi**2))
#         term *= np.sin((2*n+1) * np.pi * x / 2) * np.exp(-((2*n+1)**2 * np.pi**2 / 4) * t)
#         u_analytical += term

#     return u_analytical

# # Plot Analytical Solution at a specific time (e.g., T/2)
# t_analytical = T / 2
# u_analytical = analytical_solution(x_values, t_analytical, N_terms)

# # Plot the result
# plt.figure(figsize=(8, 6))
# plt.plot(x_values, u_analytical, label='Analytical Solution')
# plt.title('Analytical Solution at t = {:.3f}'.format(t_analytical))
# plt.xlabel('Position (x)')
# plt.ylabel('Temperature (u)')
# plt.legend()
# plt.show()






import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0  # Length of the rod
T = 0.5  # Total simulation time
Nx = 100  # Number of spatial grid points
Nt = 500  # Number of time steps
N_terms = 20  # Number of terms in the Fourier series

# Spatial and temporal discretization
dx = L / (Nx - 1)
dt = T / Nt

# Initial condition
def initial_condition(x):
    return x**3 - 3*x

# Initialize solution matrix
u = np.zeros((Nx, Nt+1))

# Set initial condition
u[:, 0] = initial_condition(np.linspace(0, L, Nx))

# Finite difference coefficients
A = -1 / dx**2
B = 1 / dt + 2 / dx**2
C = -1 / dx**2

# Time-stepping loop using implicit method and Thomas algorithm
for n in range(0, Nt):
    # Prepare the tridiagonal system
    D = u[:, n] / dt
    D[0] = 0  # Boundary condition
    D[-1] = 0  # Boundary condition

    # Thomas algorithm
    alpha = np.zeros(Nx-2)
    beta = np.zeros(Nx-2)

    for i in range(1, Nx-1):
        denominator = B + C * alpha[i-1]
        alpha[i-1] = -A / denominator
        beta[i-1] = (D[i] - C * beta[i-1]) / denominator

    # Backward substitution
    for i in range(Nx-3, -1, -1):
        u[i+1, n+1] = alpha[i] * u[i+2, n+1] + beta[i]

# Plot the numerical solution at a specific time (e.g., T/2)
t_numerical = T / 2
u_numerical = u[:, int(t_numerical / T * Nt)]

# Plot the result
plt.figure(figsize=(8, 6))
plt.plot(np.linspace(0, L, Nx), u_numerical, label='Numerical Solution')
plt.title('Numerical Solution at t = {:.3f}'.format(t_numerical))
plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')
plt.legend()
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

# # Analytical Solution
# u_analytical = np.zeros((Nx, Nt+1))
# for j, t in enumerate(t_values): 
#     u_analytical[:, j] = initial_condition(x_values)  # Using the initial condition for simplicity

# ax2 = fig.add_subplot(122, projection='3d')
# ax2.plot_surface(X, T, u_analytical.T, cmap='viridis')
# ax2.set_title('Analytical Solution')
# ax2.set_xlabel('Position (x)')
# ax2.set_ylabel('Time (t)')
# ax2.set_zlabel('Temperature (u)')

# plt.show()
