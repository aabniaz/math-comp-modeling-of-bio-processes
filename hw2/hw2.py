import numpy as np
import matplotlib.pyplot as plt

Nx = 100  
Nt = 500
L = 1.0
T = 0.5
deltax = L / (Nx - 1)
deltat = T / Nt

def heat_conduction_solver(Nx, Nt, deltax, deltat):
    A = -1 / (deltax**2)
    B = 1 / deltat + 2 / (deltax**2)
    C = -1 / (deltax**2)

    alpha = np.zeros(Nx)
    beta = np.zeros(Nx)

    alpha[1] = 0
    beta[1] = 0

    # u^n_i (initial condition)
    u_n = np.linspace(0, 1, Nx)**3 - 3 * np.linspace(0, 1, Nx)

    for i in range(1, Nx-1):
        D = u_n[i] / deltat
        denominator = B + C * alpha[i]
        alpha[i+1] = -A / denominator
        beta[i+1] = (D - C * beta[i]) / denominator

    U_N = beta[Nx-1] / (1 - alpha[Nx-1])

    # U at the last point
    U = np.zeros(Nx)
    U[Nx-1] = U_N

    # U for interior points
    for i in range(Nx-1, 1, -1):
        U[i-1] = alpha[i] * U[i] + beta[i]

    u_analytical = np.zeros(Nx)
    for n in range(1000):  
        C = 2 * (-1)**n * (4 / ((2*n+1) * np.pi) + 48 / ((2*n+1)**3 * np.pi**3) + 96 / ((2*n+1)**4 * np.pi**4))
        X = np.sin((2*n+1) * np.pi * np.linspace(0, 1, Nx) / 2)
        T_analytical = np.exp(-((2*n+1)**2 * np.pi**2 * deltat * Nt) / 4)
        u_analytical += -C * X * T_analytical

    max_difference = np.max(np.abs(U - u_analytical))

    x = np.linspace(0, 1, Nx)
    print(f'Maximum Difference: {max_difference}')

    plt.figure(figsize=(8, 6))
    plt.plot(x, U, label='Numerical Solution')
    plt.plot(x, u_analytical, label='Analytical Solution', linestyle='dashed')
    plt.xlabel('x')
    plt.ylabel('U(x, T)')
    plt.legend()
    plt.show()

heat_conduction_solver(Nx, Nt, deltax, deltat)