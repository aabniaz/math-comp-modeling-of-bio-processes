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

        if alpha**2 * deltat / deltax**2 > 0.5:
            print("Stability condition is not satisfied")

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

U_analytical = analytical(x, t=time_e, N_terms=N_terms)

print(f"Maximum difference (explicit): {np.max(np.abs(U_analytical - U_explicit))}")

plt.figure(figsize=(12, 8))
plt.plot(x, U_explicit, label='Explicit', linestyle='--')
# plt.plot(x, U_analytical, label='Analytical Solution', linestyle='dashed')
plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')
plt.legend()
plt.show()
