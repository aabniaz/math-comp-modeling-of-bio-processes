import numpy as np
import matplotlib.pyplot as plt

def analytical(x, t, N): 
    U_analytical = 0
    for n in range(1, N + 1): 
        C = 2 * (-1) ** n * (4 / ((2 * n - 1) * np.pi) + 48 / ((2 * n - 1) ** 3 * np.pi ** 3) + 96 / ( 
                    (2 * n - 1)  **4 * np.pi ** 4)) 
        X = np.sin((2 * n - 1) * np.pi * x / 2) 
        T_analytical = np.exp(-((2 * n - 1) ** 2 * np.pi ** 2 * t) / 4) 
        U_analytical += C * X * T_analytical 
    return U_analytical

iter = 0
n = 11
T = 0.1
L = 1
a = 1
dt = 0.0001
dx = 0.1
eps = 1e-6

x = np.linspace(0, L, n)
t = np.linspace(0, T, int(T / dt) + 1)

m = len(t)
n = len(x)
U = x * (x**2 - 3)

A = -a / dx**2
B = 1 / dt + 2 * a / dx**2
C = -a / dx**2
D = np.zeros(n)

alpha = np.zeros(n)
beta = np.zeros(n)
alpha[0] = 0  
beta[0] = 0

#implicit
while True:
    D = U / dt

    for i in range(1, n-1):
        alpha[i+1] = -A / (B + C*alpha[i])
        beta[i+1] = (D[i] - C*beta[i]) / (B + C*alpha[i])

    U_new = np.zeros(n)
    U_new[0] = 0
    U_new[-1] = U_new[-2]

    for i in range(n-2, 0, -1):
        U_new[i] = alpha[i+1] * U_new[i+1] + beta[i+1]

    max_diff = np.max(np.abs(U_new - U))

    if max_diff < eps:
        break

    U[:] = U_new
    iter += 1

# explicit method
new_T = iter * dt
x1 = np.arange(0, L + dx, dx)
t1 = np.arange(0, new_T + dt, dt)
U1 = np.zeros((len(t1), len(x1)))

U1[0, :] = x1 * (x1**2 - 3)
U1[:, 0] = 0
U1[:, -1] = U1[:, -2]

for i in range(len(t1)-1):
    for j in range(1, len(x1)-1):
        U1[i+1, j] = U1[i, j] + a * dt / dx**2 * (U1[i, j+1] - 2 * U1[i, j] + U1[i, j-1])


x_analytical = np.linspace(0, L, n)
u_analytical = analytical(x_analytical, new_T, iter)

print("Number of iterations:", iter)
print("Thomas:\n", U_new)
print("Euler:\n", U1[-1, :])
print("Analytical:\n", u_analytical)

x_val = np.linspace(0, L, n)
plt.plot(x_val, U_new, label='Numerical Thomas') 
plt.plot(x1, U1[-1, :], label='Numerical Euler')
# plt.plot(x_analytical, u_analytical, label='Analytical')
plt.legend()
plt.show()

print('Max difference =', np.max(np.abs(U_new - U1)))
# print('Max error for Euler =', np.max(np.abs(U1[-1, :] - u_analytical)))

print(a * dt / dx**2 <= 0.5)
print(a * dt / dx**2)
