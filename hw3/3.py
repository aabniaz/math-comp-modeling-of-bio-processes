# import numpy as np
# import matplotlib.pyplot as plt

# def implicit(Nx, Nt, deltax, deltat):
#     A = -1 / (deltax**2)
#     B = 1 / deltat + 2 / (deltax**2)
#     C = -1 / (deltax**2)

#     alpha = np.zeros(Nx)
#     beta = np.zeros(Nx)

#     alpha[1] = 0
#     beta[1] = 0

#     # u^n_i (initial condition)
#     u_n = np.linspace(0, 1, Nx)**3 - 3 * np.linspace(0, 1, Nx)

#     U_implicit = np.zeros((Nt, Nx))
#     U_implicit[0, :] = u_n
    
#     for i in range(1, Nx-1):
#         D = u_n[i] / deltat
#         denominator = B + C * alpha[i]
#         alpha[i+1] = -A / denominator
#         beta[i+1] = (D - C * beta[i]) / denominator

#     U_implicit[-1, :] = beta[Nx-1] / (1 - alpha[Nx-1])
#     U = np.zeros(Nx)
#     U[Nx-1] = U_implicit[-1, Nx-1]

#     # U for interior points
#     for i in range(Nx-1, 1, -1):
#         U[i-1] = alpha[i] * U[i] + beta[i]

#     return U

# def analytical(x, t, N_terms):
#     U_analytical = np.zeros_like(x)
#     for n in range(1, N_terms + 1):
#         C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
#                     (2 * n - 1) ** 4 * np.pi ** 4))
#         X = np.sin((2 * n - 1) * np.pi * x / 2)
#         T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
#         U_analytical += C * X * T_analytical
#     return U_analytical

# # def explicit(Nx, Nt, deltax, deltat, alpha):
# #     U_explicit = np.zeros((Nt, Nx))
# #     x_values = np.linspace(0, 1, Nx)
# #     U_explicit[0, :] = x_values**3 - 3 * x_values

# #     for n in range(Nt - 1):
# #         for i in range(1, Nx - 1):
# #             U_explicit[n + 1, i] = U_explicit[n, i] + (alpha**2 * deltat / deltax**2) * (
# #                 U_explicit[n, i+1] - 2 * U_explicit[n, i] + U_explicit[n, i-1]
# #             )

# #         # boundary conditions after the cycle by x
# #         U_explicit[n + 1, 0] = 0  # U(x=0, t) = 0
# #         U_explicit[n + 1, -1] = U_explicit[n, -2]  # U_x(x=1, t) = 0 

# #     if alpha**2 * deltat / deltax**2 > 0.5:
# #         print("Stability condition is not satisfied!")

# #     return U_explicit

# Nx = 100  
# Nt = 500
# alpha = 1
# deltax = 0.1
# deltat = 0.0001
# epsilon = 1e-6  
# N_terms = 20

# U_implicit = implicit(Nx, Nt, deltax, deltat)
# U_explicit = explicit(Nx, Nt, deltax, deltat, alpha)
# x_values = np.linspace(0, 1, Nx)
# t_values = np.linspace(0, 0.5, Nt)
# U_analytical = np.zeros((Nx, Nt))
# for j, t in enumerate(t_values):
#     U_analytical[:, j] = analytical(x_values, t, N_terms)

# plt.figure(figsize=(12, 8))
# # plt.plot(x_values, U_implicit, label='Implicit', linestyle='--')
# plt.plot(x_values, U_explicit[-1, :], label='Explicit', linestyle='--')
# plt.plot(x_values, U_analytical[:, -1], label='Analytical', linestyle='-')
# plt.xlabel('Position (x)')
# plt.ylabel('Temperature (u)')
# plt.legend()
# plt.show()
# for j, t in enumerate(t_values):
#     U_analytical[:, j] = analytical(x_values, t, N_terms)
# if U_explicit.shape[0] == U_analytical.shape[1]:
#     last_time_moment = -1
#     max_difference = np.max(np.abs(U_explicit[last_time_moment, :] - U_analytical[:, last_time_moment]))
#     print("Max Difference at the Last Time Moment:", max_difference)
# else:
#     print("Shapes mismatch. Unable to compute max difference.")
# print(np.max(np.abs(U_explicit[last_time_moment, :] - U_analytical[:, last_time_moment])) <= deltat + deltax**2)
# print(deltat + deltax**2)






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
    U_analytical = np.zeros_like(x)
    for n in range(1, N_terms + 1):
        C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
                    (2 * n - 1) ** 4 * np.pi ** 4))
        X = np.sin((2 * n - 1) * np.pi * x / 2)
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
        U_analytical += C * X * T_analytical
    return U_analytical

# def explicit(U, deltax, deltat, alpha=1, eps=1e-6):
#     N = len(U)
#     U_init = np.copy(U)
#     U_new = np.zeros(N)

#     iter = 0
#     while True:
#         U_new[1:N-1] = U_init[1:N-1] + alpha**2 * deltat / deltax**2 * (U_init[2:N] - 2*U_init[1:N-1] + U_init[0:N-2])
#         U_new[0] = 0
#         U_new[N-1] = U_new[N-2]

#         maximum = np.max(np.abs(U_new - U_init))
#         U_init[:] = U_new[:]
#         iter += 1

#         if maximum <= eps:
#             break

#     return U_new, iter*deltat, iter

def explicit(U, deltax, deltat, alpha=1, eps=1e-6):
    N = len(U)
    U_init = np.copy(U)
    U_new = np.zeros(N)

    iter = 0
    while True:
        U_new[1:N-1] = U_init[1:N-1] + alpha**2 * deltat / deltax**2 * (U_init[2:N] - 2*U_init[1:N-1] + U_init[0:N-2])
        U_new[0] = 0
        U_new[N-1] = U_new[N-2]

        maximum = np.max(np.abs(U_new - U_init))
        U_init[:] = U_new[:]
        iter += 1

        if maximum <= eps:
            break

    return U_new, iter*deltat, iter

# Nx = 100  
# Nt = 500
# deltax = 0.1
# deltat = 0.0001
# alpha = 1

# # Initial conditions
# U_init = np.linspace(0, 1, Nx)**3 - 3 * np.linspace(0, 1, Nx)

# # Solving implicit and explicit schemes
# U_implicit = implicit(Nx, Nt, deltax, deltat)
# U_explicit, _, _ = explicit(U_init, deltax, deltat, alpha)

# # Analytical solution
# x_values = np.linspace(0, 1, Nx)
# t_values = np.linspace(0, 0.5, Nt)
# U_analytical = np.zeros((Nx, Nt))
# for j, t in enumerate(t_values):
#     U_analytical[:, j] = analytical(x_values, t, N_terms=20)

# # Plotting
# plt.figure(figsize=(12, 8))
# plt.plot(x_values, U_implicit, label='Implicit', linestyle='--')
# plt.plot(x_values, U_explicit, label='Explicit', linestyle='--')
# plt.plot(x_values, U_analytical[:, -1], label='Analytical', linestyle='-')
# plt.xlabel('Position (x)')
# plt.ylabel('Temperature (u)')
# plt.legend()
# plt.show()

# # Calculate and print maximum difference
# last_time_moment = -1
# max_difference = np.max(np.abs(U_explicit - U_analytical[:, last_time_moment]))
# print("Max Difference at the Last Time Moment:", max_difference)
# print(np.max(np.abs(U_explicit - U_analytical[:, last_time_moment])) <= deltat + deltax**2)


# Parameters for implicit scheme
Nx_implicit = 100  
Nt_implicit = 500
deltax_implicit = 0.1
deltat_implicit = 0.0001

# Parameters for explicit scheme
Nx_explicit = 100  
Nt_explicit = 500
deltax_explicit = 0.1
deltat_explicit = 0.0001
alpha_explicit = 1

# Initial conditions
U_init_explicit = np.linspace(0, 1, Nx_explicit)**3 - 3 * np.linspace(0, 1, Nx_explicit)

# Solving explicit scheme
U_explicit, time_explicit, iter_explicit = explicit(U_init_explicit, deltax_explicit, deltat_explicit, alpha_explicit)

# Analytical solution
x_values_explicit = np.linspace(0, 1, Nx_explicit)
t_values_explicit = np.linspace(0, 0.5, Nt_explicit)
U_analytical_explicit = np.zeros((Nx_explicit, Nt_explicit))
for j, t in enumerate(t_values_explicit):
    U_analytical_explicit[:, j] = analytical(x_values_explicit, t, N_terms=20)

# Print maximum difference and check condition
max_difference_explicit = np.max(np.abs(U_analytical_explicit[:, -1] - U_explicit))
print(f"Maximum difference: {max_difference_explicit}")
print(max_difference_explicit <= deltat_explicit + deltax_explicit**2)


# Plotting
plt.figure(figsize=(12, 8))
plt.plot(x_values_explicit, U_explicit, label='Explicit', linestyle='--')
plt.plot(x_values_explicit, U_analytical_explicit[:, -1], label='Analytical', linestyle='-')
plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')
plt.legend()
plt.show()