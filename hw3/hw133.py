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

def Thomas_algorithm(u_n, deltax, deltat, Nx):
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

dx = 0.1
dt = 0.001
eps = 1e-6
start_x, end_x = (0, 1)
Nx = int((end_x - start_x) / dx) + 1
x = start_x + np.arange(start=0, stop=Nx) * dx
U_init = np.zeros_like(x)
U_init[:] = x**3-3*x

# Checking stability condition for Simple Iterative Method
if 1 * dt / dx**2 <= 0.5:
    print("Simple Iterative Method is stable")
else:
    print("Simple Iterative Method is not stable")

# Time-stepping loop
for t_step in range(1, 1001):
    t = t_step * dt
    
    U_analytical = analytical(x, t, N=100)
    U_init = Thomas_algorithm(U_init, dx, dt, Nx)
    U_explicit, time_s, iter_s = explicit(U_init, dx, dt)
    exp_analytical = analytical(x, t=time_s, N=100)
    
    if t_step % 100 == 0:
        plt.figure(figsize=(12, 6))        
        plt.subplot(1, 2, 1)
        plt.grid()
        plt.plot(x, U_init, label=f"Implicit at t={t}")
        plt.plot(x, U_analytical, ls="--", label=f"Analytical at t={t}")
        plt.title("Heat equation - Implicit Method")
        plt.xlabel("x")
        plt.ylabel("U(x)")
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.grid()
        plt.plot(x, U_explicit, label=f"Simple method at t={time_s}")
        plt.plot(x, exp_analytical, ls="--", label=f"Analytical for (Simple) at t={time_s}")
        plt.title("Heat equation - Simple Iterative Method")
        plt.xlabel("x")
        plt.ylabel("U(x)")
        plt.legend()

        plt.tight_layout()
        plt.show()

max_diff_thomas = np.max(np.abs(U_analytical - U_init))
max_diff_simple = np.max(np.abs(exp_analytical - U_explicit))

print("Maximum difference (Thomas algorithm):", max_diff_thomas)
print("Maximum difference (Simple Iterative Method):", max_diff_simple)



# import numpy as np
# import matplotlib.pyplot as plt
# from numpy import pi, exp, sin

# def Analytical_solution(x, t, m=100):
#     λ = lambda n: (2*n - 1)*pi/2
#     A_n = lambda n: -12*((-1)**n + 1/λ(n)) / λ(n)**3
#     A_n = lambda n: 2*(-1)**n * (2 / λ(n)  + 24 / λ(n)**3 + 48 / λ(n)**4)

#     F = np.zeros_like(x)
    
#     for n in range(1, m+1):
#         F += A_n(n) * exp(-λ(n)**2*t) * sin(λ(n)*x)
    
#     return F

# dx = 0.1
# dt = 0.001
# eps = 1e-6
# a2 = 1
# stop_iteration = 3e5
# start_x, end_x = (0, 1)
# N = int((end_x - start_x) / dx) + 1
# x = start_x + np.arange(start=0, stop=N) * dx
# U_old = np.zeros_like(x)
# U_old[:] =x**3-3*x

# if a2 * dt / dx**2 <= 0.5:
#     print("It is stable")
# else:
#     print("It is not stable")

# def Thomas_algorithm(U, dx, dt, a1, b1, eps=1e-6, stop_iteration=1e5):
#     N = len(U)
#     U_old = np.copy(U)
#     U_new = np.zeros_like(U)
#     alpha = np.zeros_like(U)
#     beta = np.zeros_like(U)

#     A = -1/dx**2
#     B = 1/dt + 2/dx**2
#     C = -1/dx**2
#     D = np.zeros_like(U)

#     comparison = np.abs(B) >= np.abs(A) + np.abs(C)
#     print("Thomas algorithm")
    
#     if comparison.all():
#         print("it is stable")
#     else:
#         print("It is not stable")
#         return

#     iteration = 0
#     maximum = 1
    
#     while maximum > eps and iteration < stop_iteration:
#         D[:] = U_old/dt
#         alpha[1] = a1
#         beta[1] = b1
        
#         for i in range(1, N-1):
#             alpha[i+1] = -A/(B + C*alpha[i])
#             beta[i+1] = (D[i] - C*beta[i])/(B + C*alpha[i])

#         U_new[N-1] = 0
        
#         for i in range(N-2, 0, -1):
#             U_new[i] = alpha[i+1]*U_new[i+1] + beta[i+1]
        
#         U_new[0] = U_new[1]
#         maximum = np.max(np.abs(U_new - U_old))
#         U_old[:] = U_new[:]
#         iteration += 1

#     return U_new, iteration*dt, iteration

# def Simple_Iterative_Method(U, dx, dt, a2=1, eps=1e-6, stop_iteration=1e5):
#     N = len(U)
#     U_old = np.copy(U)
#     U_new = np.zeros(N)
    
#     iteration = 0
#     maximum = 1
    
#     while maximum > eps and iteration < stop_iteration:
#         U_new[1:N-1] = U_old[1:N-1] + a2*dt/dx**2 * (U_old[2:N] - 2*U_old[1:N-1] + U_old[0:N-2])
#         U_new[0] = 0
#         U_new[N-1] = U_new[N-2]
        
#         maximum = np.max(np.abs(U_new - U_old))
#         U_old[:] = U_new[:]
#         iteration += 1

#     return U_new, iteration*dt, iteration

# TM, time_t, iter_t = Thomas_algorithm(U_old, dx, dt, a1=0, b1=0)
# U_analytical = Analytical_solution(x, t=time_t)
# U_explicit, time_s, iter_s = Simple_Iterative_Method(U_old, dx, dt)
# exp_analytical = Analytical_solution(x, t=time_s)

# print(time_t, time_s)

# plt.grid()
# plt.plot(x, TM, label=f"Thomas algorithm at {time_t}")
# plt.plot(x, U_analytical, ls="--", label=f"Analytical for (Thomas) at {time_t}")
# plt.plot(x, U_explicit, label=f"Simple method at {time_s}")
# plt.plot(x, exp_analytical, ls="--", label=f"Analytical for (Simple) at {time_s}")
# plt.title("Heat equation")
# plt.xlabel("x")
# plt.ylabel("U(x)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# print("Thomas algorithm")
# print(f"Maximum difference: {np.max(np.abs(U_analytical - TM)):.9f}")
# print(f"Number of iterations: {iter_t}\n")

# print("Simple iterative method")
# print(f"Maximum difference: {np.max(np.abs(exp_analytical - U_explicit)):.9f}")
# print(f"Number of iterations: {iter_s}\n")

# print("Difference between two methods")
# print(f"Maximum difference: {np.max(np.abs(TM - U_explicit)):.9f}")
# print(np.max(np.abs(exp_analytical - U_explicit)) <= dt + dx**2)



# import numpy as np
# import matplotlib.pyplot as plt
# from numpy import pi, exp, sin

# def Analytical_solution(x, t, m=100):
#     λ = lambda n: (2*n - 1)*pi/2
#     A_n = lambda n: 2*(-1)**n * (2 / λ(n)  + 24 / λ(n)**3 + 48 / λ(n)**4)

#     F = np.zeros_like(x)
    
#     for n in range(1, m+1):
#         F += A_n(n) * exp(-λ(n)**2*t) * sin(λ(n)*x)
    
#     return F

# def Thomas_algorithm(U, dx, dt, a1, b1, eps=1e-6, stop_iteration=1e5):
#     N = len(U)
#     U_old = np.copy(U)
#     U_new = np.zeros_like(U)
#     alpha = np.zeros_like(U)
#     beta = np.zeros_like(U)

#     A = -1/dx**2
#     B = 1/dt + 2/dx**2
#     C = -1/dx**2
#     D = np.zeros_like(U)

#     comparison = np.abs(B) >= np.abs(A) + np.abs(C)
#     print("Thomas algorithm")
    
#     if comparison.all():
#         print("it is stable")
#     else:
#         print("It is not stable")
#         return

#     iteration = 0
#     maximum = 1
    
#     while maximum > eps and iteration < stop_iteration:
#         D[:] = U_old/dt
#         alpha[1] = a1
#         beta[1] = b1
        
#         for i in range(1, N-1):
#             alpha[i+1] = -A/(B + C*alpha[i])
#             beta[i+1] = (D[i] - C*beta[i])/(B + C*alpha[i])

#         U_new[N-1] = 0
        
#         for i in range(N-2, 0, -1):
#             U_new[i] = alpha[i+1]*U_new[i+1] + beta[i+1]
        
#         U_new[0] = U_new[1]
#         maximum = np.max(np.abs(U_new - U_old))
#         U_old[:] = U_new[:]
#         iteration += 1

#     return U_new, iteration*dt, iteration

# def Simple_Iterative_Method(U, dx, dt, a2=1, eps=1e-6, stop_iteration=1e5):
#     N = len(U)
#     U_old = np.copy(U)
#     U_new = np.zeros(N)
    
#     iteration = 0
#     maximum = 1
    
#     while maximum > eps and iteration < stop_iteration:
#         U_new[1:N-1] = U_old[1:N-1] + a2*dt/dx**2 * (U_old[2:N] - 2*U_old[1:N-1] + U_old[0:N-2])
#         U_new[0] = 0
#         U_new[N-1] = U_new[N-2]
        
#         maximum = np.max(np.abs(U_new - U_old))
#         U_old[:] = U_new[:]
#         iteration += 1

#     return U_new, iteration*dt, iteration

# dx = 0.1
# dt = 0.001
# eps = 1e-6
# a2 = 1
# stop_iteration = 3e5
# start_x, end_x = (0, 1)
# N = int((end_x - start_x) / dx) + 1
# x = start_x + np.arange(start=0, stop=N) * dx
# U_old = np.zeros_like(x)
# U_old[:] =x**3-3*x

# if a2 * dt / dx**2 <= 0.5:
#     print("It is stable")
# else:
#     print("It is not stable")

# TM, time_t, iter_t = Thomas_algorithm(U_old, dx, dt, a1=0, b1=0)
# U_analytical = Analytical_solution(x, t=time_t)
# U_explicit, time_s, iter_s = Simple_Iterative_Method(U_old, dx, dt)
# exp_analytical = Analytical_solution(x, t=time_s)

# print(time_t, time_s)

# # Subplots for Implicit and Analytical
# plt.figure(figsize=(12, 6))
# plt.subplot(1, 2, 1)
# plt.grid()
# plt.plot(x, TM, label=f"Thomas algorithm at {time_t}")
# plt.plot(x, U_analytical, ls="--", label=f"Analytical for (Thomas) at {time_t}")
# plt.title("Heat equation - Implicit Method")
# plt.xlabel("x")
# plt.ylabel("U(x)")
# plt.legend()

# # Subplots for Explicit and Analytical
# plt.subplot(1, 2, 2)
# plt.grid()
# plt.plot(x, U_explicit, label=f"Simple method at {time_s}")
# plt.plot(x, exp_analytical, ls="--", label=f"Analytical for (Simple) at {time_s}")
# plt.title("Heat equation - Simple Iterative Method")
# plt.xlabel("x")
# plt.ylabel("U(x)")
# plt.legend()

# plt.tight_layout()
# plt.show()

# print("Thomas algorithm")
# print(f"Maximum difference: {np.max(np.abs(U_analytical - TM)):.9f}")
# print(f"Number of iterations: {iter_t}\n")

# print("Simple iterative method")
# print(f"Maximum difference: {np.max(np.abs(exp_analytical - U_explicit)):.9f}")
# print(f"Number of iterations: {iter_s}\n")

# print("Difference between two methods")
# print(f"Maximum difference: {np.max(np.abs(TM - U_explicit)):.9f}")
# print(np.max(np.abs(exp_analytical - U_explicit)) <= dt + dx**2)