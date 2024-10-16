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
    
    for j in range(1, Nt):
        D = U_implicit[j-1, 1:Nx-1] / deltat
        denominator = B + C * alpha[1:Nx-1]
        alpha[2:Nx] = -A / denominator
        beta[2:Nx] = (D - C * beta[1:Nx-1]) / denominator

        U_implicit[j, -1] = beta[Nx-1] / (1 - alpha[Nx-1])
        U_implicit[j, 1:Nx-1] = alpha[2:Nx] * U_implicit[j, 1:Nx-1] + beta[2:Nx]

    return U_implicit


def analytical(x, t, N_terms):
    U_analytical = np.zeros_like(x)
    for n in range(1, N_terms + 1):
        C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
                    (2 * n - 1) ** 4 * np.pi ** 4))
        X = np.sin((2 * n - 1) * np.pi * x / 2)
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
        U_analytical += C * X * T_analytical
    return U_analytical

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
    
Nx = 50  
Nt = 500
deltax = 0.1
deltat = 0.001
eps = 1e-6
alpha = 1
stop_it = 5000
start_x, end_x = (0, 1)
N = 50
x = np.linspace(start_x, end_x, N)
U_init = x**3 - 3*x 
N_terms = 20

U_implicit = implicit(Nx, Nt, deltax, deltat)
U_explicit, time_e, iter_e = explicit(U_init, deltax, deltat, alpha, eps)

# Collecting solutions for each time step for explicit method
# U_explicit_solution = np.zeros((Nt, Nx))
# U_explicit_solution[0, :] = U_init
# for i in range(1, Nt):
#     U_explicit_solution[i, :], _, _ = explicit(U_explicit_solution[i-1, :], deltax, deltat, alpha, eps)



U_analytical = analytical(x, t=time_e, N_terms=20)

print(f"Maximum difference: {np.max(np.abs(U_analytical - U_explicit))}")
print(np.max(np.abs(U_analytical - U_explicit)) <= deltat + deltax**2)

print("U_explicit:", U_explicit)
print("U_analytical:", U_analytical)

x_values = np.linspace(0, 1, Nx)
t_values = np.linspace(0, 0.5, Nt)
U_analytical = np.zeros((Nx, Nt))
for j, t in enumerate(t_values):
    U_analytical[:, j] = analytical(x_values, t, N_terms)

plt.figure(figsize=(12, 8))

# Find the index corresponding to t = 0.05
index_t_50 = int(0.09 / deltat)

# Plot Implicit at t = 0.05
plt.plot(x_values, U_implicit[index_t_50, :], label='Implicit (t = 0.05)', linestyle='--')

# Plot Analytical
U_analytical_final = analytical(x_values, t=time_e, N_terms=20)
plt.plot(x_values, U_analytical_final, label='Analytical', linestyle='-')

# Plot Explicit
plt.plot(x_values, U_explicit, label='Explicit', linestyle='--')

plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')



plt.legend()
plt.show()


# import matplotlib.pyplot as plt
# import numpy as np

# # Given data
# U_explicit = np.array([0.0, -0.00030872, -0.00061711, -0.00092486, -0.00123164, -0.00153712,
#                       -0.00184099, -0.00214294, -0.00244263, -0.00273976, -0.00303402, -0.0033251,
#                       -0.00361269, -0.00389649, -0.0041762, -0.00445153, -0.0047222, -0.00498791,
#                       -0.00524839, -0.00550336, -0.00575256, -0.00599573, -0.00623261, -0.00646295,
#                       -0.00668652, -0.00690307, -0.00711238, -0.00731423, -0.00750841, -0.00769471,
#                       -0.00787294, -0.00804292, -0.00820446, -0.00835739, -0.00850156, -0.00863681,
#                       -0.008763, -0.00888, -0.00898769, -0.00908595, -0.00917468, -0.00925379,
#                       -0.00932319, -0.00938281, -0.00943259, -0.00947248, -0.00950243, -0.00952242,
#                       -0.00953242, -0.00953242])

# U_analytical = np.array([0.0, -9.48343567e-56, -1.89571265e-55, -2.84113375e-55, -3.78363541e-55,
#                          -4.72224911e-55, -5.65601039e-55, -6.58395974e-55, -7.50514362e-55,
#                          -8.41861546e-55, -9.32343660e-55, -1.02186773e-54, -1.11034176e-54,
#                          -1.19767484e-54, -1.28377722e-54, -1.36856044e-54, -1.45193737e-54,
#                          -1.53382234e-54, -1.61413120e-54, -1.69278143e-54, -1.76969222e-54,
#                          -1.84478452e-54, -1.91798119e-54, -1.98920700e-54, -2.05838877e-54,
#                          -2.12545540e-54, -2.19033798e-54, -2.25296984e-54, -2.31328663e-54,
#                          -2.37122635e-54, -2.42672948e-54, -2.47973898e-54, -2.53020038e-54,
#                          -2.57806183e-54, -2.62327415e-54, -2.66579088e-54, -2.70556832e-54,
#                          -2.74256561e-54, -2.77674473e-54, -2.80807056e-54, -2.83651090e-54,
#                          -2.86203654e-54, -2.88462124e-54, -2.90424180e-54, -2.92087805e-54,
#                          -2.93451290e-54, -2.94513235e-54, -2.95272547e-54, -2.95728446e-54,
#                          -2.95880465e-54])

# # Plotting absolute values
# plt.plot(-np.abs(U_explicit), label='|U_explicit|')
# plt.plot(U_analytical, label='U_analytical')
# plt.xlabel('Index')
# plt.ylabel('Absolute Values')
# plt.title('Comparison of |U_explicit| and |U_analytical|')
# plt.legend()
# plt.show()



# U_explicit = np.array([0.0, -0.00030872, -0.00061711, -0.00092486, -0.00123164, -0.00153712,
#                       -0.00184099, -0.00214294, -0.00244263, -0.00273976, -0.00303402, -0.0033251,
#                       -0.00361269, -0.00389649, -0.0041762, -0.00445153, -0.0047222, -0.00498791,
#                       -0.00524839, -0.00550336, -0.00575256, -0.00599573, -0.00623261, -0.00646295,
#                       -0.00668652, -0.00690307, -0.00711238, -0.00731423, -0.00750841, -0.00769471,
#                       -0.00787294, -0.00804292, -0.00820446, -0.00835739, -0.00850156, -0.00863681,
#                       -0.008763, -0.00888, -0.00898769, -0.00908595, -0.00917468, -0.00925379,
#                       -0.00932319, -0.00938281, -0.00943259, -0.00947248, -0.00950243, -0.00952242,
#                       -0.00953242, -0.00953242])

# plt.plot(np.abs(U_explicit), label='|U_explicit|')
# plt.xlabel('Index')
# plt.ylabel('Absolute Values')
# plt.title('Comparison of |U_explicit| and |U_analytical|')
# plt.legend()
# plt.show()
