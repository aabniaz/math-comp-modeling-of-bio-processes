import numpy as np
import matplotlib.pyplot as plt

def analytical(x, t, N=100):
    U_analytical = np.zeros_like(x)
    for n in range(1, N + 1):
        C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / (
                    (2 * n - 1) ** 4 * np.pi ** 4))
        X = np.sin((2 * n - 1) * np.pi * x / 2)
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4)
        U_analytical += C * X * T_analytical
    return U_analytical

def implicit(u_n, deltax, deltat, Nx):
    A = -1 / (deltax**2)
    B = 1 / deltat + 2 / (deltax**2)
    C = -1 / (deltax**2)

    alpha = np.zeros(Nx)
    beta = np.zeros(Nx)

    alpha[1] = 0
    beta[1] = 0

    # u^n_i (initial condition)
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

    return U

def explicit(U, deltax, deltat, eps=1e-6):
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

dx = 0.1
dt = 0.001
eps = 1e-6
alpha = 1
start_x, end_x = (0, 1)
Nx = int((end_x - start_x) / dx) + 1
x = start_x + np.arange(start=0, stop=Nx) * dx
U_init = np.zeros_like(x)
U_init[:] = x**3-3*x

if alpha**2 * dt / dx**2 <= 0.5:
    print("stability condition is satisfied")
else:
    print("stability condition is not satisfied")

for t_step in range(1, 1001):
    t = t_step * dt
    
    U_analytical = analytical(x, t, N=100)
    U_init = implicit(U_init, dx, dt, Nx)
    U_explicit, time_s, iter_s = explicit(U_init, dx, dt)
    exp_analytical = analytical(x, t=time_s, N=100)
    
    if t_step % 100 == 0:
        plt.figure(figsize=(12, 6))        
        plt.subplot(1, 2, 1)
        plt.grid()
        plt.plot(x, U_init, label=f"implicit at t={t}")
        plt.plot(x, U_analytical, ls="--", label=f"analytical at t={t}")
        plt.title("heat equation - implicit method")
        plt.xlabel("x")
        plt.ylabel("U(x)")
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.grid()
        plt.plot(x, U_explicit, label=f"explicit at t={time_s}")
        plt.plot(x, exp_analytical, ls="--", label=f"analytical for (explicit) at t={time_s}")
        plt.title("heat equation - explicit method")
        plt.xlabel("x")
        plt.ylabel("U(x)")
        plt.legend()

        plt.tight_layout()
        # plt.show()
        #plt.savefig(f'figure_{t_step}.png')  

# print(f'max difference (implicit): {np.max(np.abs(U_analytical - U_init))}')
print(f'max difference (explicit): {np.max(np.abs(exp_analytical - U_explicit))}')