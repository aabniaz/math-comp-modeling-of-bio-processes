import numpy as np
import matplotlib.pyplot as plt

start_x, end_x = (0, 1)
start_y, end_y = (0, 1)
a2 = 1
dt = 0.01
eps = 1e-6
stop_iteration = 10000
dx = 0.01
dy = 0.01
x = np.arange(0, dx+ 1, dx)
y = np.arange(0, dy + 1, dy)

M = len(x)
N = len(y)

M1 = int(0.3 * M)
M2 = int(0.7 * M)

X, Y = np.meshgrid(x, y)
# U(y, x)
U_old = np.zeros((M, N))
U_new = np.zeros((M, N))
A = np.zeros((M, N))
B = np.zeros((M, N))
C = np.zeros((M, N))
D = np.zeros((M, N))
alpha = np.zeros((M, N))
beta = np.zeros((M, N))

# Initial condition
U_old[:, :] = 0
# Boundary conditions
# left
U_old[:M1, 0] = 0 #x = 0 by y 
U_old[M1:M2, 0] = 1
U_old[M2:, 0] = 0
# right
U_old[:M1, M-1] = 1
U_old[M1:, M-1] = 0
# bottom
U_old[0, M-1] = 0 #y = 0 by x 
# top
U_old[N-1, :] = 0

iteration = 0
maximum = 1
while maximum > eps and iteration < stop_iteration:
    # Finding U^(n+1/2)
    A[0:M, 0:N] = -a2 / (2*dx**2)
    B[0:M, 0:N] = 1 / dt + a2 / dx**2
    C[0:M, 0:N] = -a2 / (2*dx**2)
    D[1:M-1, 1:N-1] = U_old[1:M-1, 1:N-1] / dt \
    + a2*(U_old[1:M-1, 2:N] - 2*U_old[1:M-1, 1:N-1] + U_old[1:M-1, 0:N-2]) \
    / (2*dx**2) \
    + a2*(U_old[2:M, 1:N-1] - 2*U_old[1:M-1, 1:N-1] + U_old[0:M-2, 1:N-1]) \
    / dy**2
    # Thomas algorithm for x
    alpha[0:M2, 1] = 0
    beta[0:M2, 1] = 0
    alpha[M2:M, 1] = 0
    beta[M2:M, 1] = 1
    # U_new[0:M, N-1] = U_old[0:M, N-1]
    for j in range(1, M-1, 1):
        for i in range(1, N-1, 1):
            alpha[j, i+1] = -A[j, i] \
            / (B[j, i] + C[j, i]*alpha[j, i])
            beta[j, i+1] = (D[j, i] - C[j, i]*beta[j, i]) \
            / (B[j, i] + C[j, i]*alpha[j, i])
        # U^(n+1/2)
        U_new[j, N-1] = U_old[j, N-1]
        for i in range(N-2, -1, -1):
            U_new[j, i] = alpha[j, i+1]*U_new[j, i+1] + beta[j, i+1]
    # Finding U^(n+1)
    A[0:M, 0:N] = -a2 / (2*dy**2)
    B[0:M, 0:N] = 1 / dt + a2 / dy**2
    C[0:M, 0:N] = -a2 / (2*dy**2)
    D[1:M-1, 1:N-1] = U_new[1:M-1, 1:N-1] / dt \
    - a2*(U_old[2:M, 1:N-1] - 2*U_old[1:M-1, 1:N-1] + U_old[0:M-2, 1:N-1]) \
        / (2*dy**2)
    # Thomas algorithm for y
    alpha[1, 0:N] = 0
    beta[1, 0:N] = 0
    for i in range(1, N-1, 1):
        for j in range(1, M-1, 1):
            alpha[j+1, i] = -A[j, i] \
            / (B[j, i] + C[j, i]*alpha[j, i])
            beta[j+1, i] = (D[j, i] - C[j, i]*beta[j, i]) \
            / (B[j, i] + C[j, i]*alpha[j, i])
        # U^(n+1)
        U_new[M-1, i] = U_old[M-1, i]
        for j in range(M-2, -1, -1):
            U_new[j, i] = alpha[j+1, i]*U_new[j+1, i] + beta[j+1, i]
    # Boundary conditions for U^(n+1)
    # left
    U_new[:M1, 0] = 0 #x = 0 by y 
    U_new[M1:M2, 0] = 1
    U_new[M2:, 0] = 0
    # right
    U_new[:M1, M-1] = 1
    U_new[M1:, M-1] = 0
    # bottom
    U_new[0, M-1] = 0 #y = 0 by x 
    # top
    U_new[N-1, :] = 0
    maximum = np.max(np.abs(U_new - U_old))
    U_old[:, :] = U_new[:, :]
    iteration += 1

print(f"Number of iterations: {iteration}")
plt.title("Fractional step method")
plt.xlabel("x")
plt.ylabel("y")
plt.contourf(X, Y, U_new)
plt.colorbar()
plt.tight_layout
plt.show()