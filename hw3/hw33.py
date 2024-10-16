import numpy as np
import matplotlib.pyplot as plt

def implicit(U, dx, dt, a1, b1, eps=1e-6, stop_it=1e5):
    N = len(U)
    U_old = np.copy(U)
    U_new = np.zeros_like(U)
    alpha = np.zeros_like(U)
    beta = np.zeros_like(U)

    A = -1/dx**2
    B = 1/dt + 2/dx**2
    C = -1/dx**2
    D = np.zeros_like(U)

    comparison = np.abs(B) >= np.abs(A) + np.abs(C)
    print("Thomas algorithm")
    if comparison.all():
        print("It is stable")
    else:
        print("It is not stable")
        return

    iteration = 0
    maximum = 1
    while maximum > eps and iteration < stop_it:
        D[:] = U_old/dt

        alpha[1] = a1
        beta[1] = b1
        for i in range(1, N-1):
            alpha[i+1] = -A/(B + C*alpha[i])
            beta[i+1] = (D[i] - C*beta[i])/(B + C*alpha[i])

        U_new[N-1] = 0
        for i in range(N-2, 0, -1):
            U_new[i] = alpha[i+1]*U_new[i+1] + beta[i+1]
        U_new[0] = U_new[1]

        maximum = np.max(np.abs(U_new - U_old))
        U_old[:] = U_new[:]
        iteration += 1

    return U_new, iteration*dt, iteration

def analytical(x, t, N_terms):
    U_analytical = np.zeros_like(x)
    for n in range(1, N_terms + 1):
        C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
                    (2 * n - 1) ** 4 * np.pi ** 4))
        X = np.sin((2 * n - 1) * np.pi * x / 2)
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
        U_analytical += C * X * T_analytical
    return U_analytical



def explicit(Nx, Nt, deltax, deltat, alpha):
    U_explicit = np.zeros((Nt, Nx))
    x_values = np.linspace(0, 1, Nx)
    U_explicit[0, :] = x_values**3 - 3 * x_values

    for n in range(Nt - 1):
        for i in range(1, Nx - 1):
            U_explicit[n + 1, i] = U_explicit[n, i] + (alpha**2 * deltat / deltax**2) * (
                U_explicit[n, i+1] - 2 * U_explicit[n, i] + U_explicit[n, i-1]
            )
            #print(f'n={n}\t i={i}\t U_explicit={U_explicit[n + 1, i]:.6f}')

            # boundary conditions after the cycle by x
            U_explicit[n + 1, 0] = 0  # U(x=0, t) = 0
            U_explicit[n + 1, -1] = U_explicit[n, -2]  # U_x(x=1, t) = 0 

    if alpha**2 * deltat / deltax**2 > 0.5:
        print("Stability condition is not satisfied!")

    # # |u^(n+1) - u^n| < epsilon
    # if np.max(np.abs(U_explicit[n + 1, :] - U_explicit[n, :])) > epsilon:
    #     print("Explicit scheme may be unstable.")

    return U_explicit

Nx = 100  
Nt = 500
alpha = 1
deltax = 0.1
deltat = 0.0001
epsilon = 1e-6  
N_terms = 20

U_implicit = implicit(Nx, Nt, deltax, deltat,b1=0)
U_explicit = explicit(Nx, Nt, deltax, deltat, alpha)
x_values = np.linspace(0, 1, Nx)
t_values = np.linspace(0, 0.5, Nt)
U_analytical = np.zeros((Nx, Nt))
for j, t in enumerate(t_values):
    U_analytical[:, j] = analytical(x_values, t, N_terms)

plt.figure(figsize=(12, 8))
plt.plot(x_values, U_implicit, label='Implicit', linestyle='--')
plt.plot(x_values, U_explicit[-1, :], label='Explicit', linestyle='--')
plt.plot(x_values, U_analytical[:, -1], label='Analytical', linestyle='-')
plt.xlabel('Position (x)')
plt.ylabel('Temperature (u)')
plt.legend()
plt.show()
#print(U_explicit.shape)
for j, t in enumerate(t_values):
    U_analytical[:, j] = analytical(x_values, t, N_terms)

if U_explicit.shape[0] == U_analytical.shape[1]:
    last_time_moment = -1
    max_difference = np.max(np.abs(U_explicit[last_time_moment, :] - U_analytical[:, last_time_moment]))
    print("Max Difference at the Last Time Moment:", max_difference)
else:
    print("Shapes mismatch. Unable to compute max difference.")
print(np.max(np.abs(U_explicit[last_time_moment, :] - U_analytical[:, last_time_moment]))<=deltat+deltax**2)
print(deltat+deltax**2)









