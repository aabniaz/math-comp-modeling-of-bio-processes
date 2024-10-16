import numpy as np
import matplotlib.pyplot as plt

a = 1
eps = 1e-6
iter = 0
max = 1
stop_iter = 10000
dx = 0.01
dy = 0.01
dt = 0.0001

x = np.arange(0, dx+1, dx)
y = np.arange(0, dy+1, dy)
M = len(x)
N = len(y)
M1 = int(0.3*M)
M2 = int(0.6*M)

A = -a/(2*dx**2)
B = 1/dt + a/dx**2
C = -a/(2*dx**2)

U_old = np.zeros((M, N))
U_new = np.zeros((M, N))

# left (x = 0 by y)
U_old[:M1, 0] = 0 
U_old[M1:M2, 0] = 1
U_old[M2:, 0] = 0
# right
U_old[:M1, M-1] = 1  
U_old[M1:, M-1] = 0
# top
U_old[N-1, :] = 0  
# bottom (y = 0 by x)
U_old[0, M-1] = 0 

while max > eps and iter < stop_iter:
    # Step 1 calculations
    D = np.zeros_like(U_old)
    alpha = np.zeros_like(U_old)
    beta = np.zeros_like(U_old)

    D[1:-1, 1:-1] = U_old[1:-1, 1:-1]/dt + a * (U_old[1:-1, 2:] - 2*U_old[1:-1, 1:-1] + U_old[1:-1, :-2])/(2*dx**2) + a*(U_old[2:, 1:-1] - 2*U_old[1:-1, 1:-1] + U_old[:-2, 1:-1])/(dy**2)

    beta[:, 1] = U_old[:, 0]

    for i in range(1, M-1):
        alpha[:, i+1] = -np.full(M, A) / (np.full(M, B) + C*alpha[:, i])
        beta[:, i+1] = (D[:, i] - C*beta[:, i]) / (np.full(M, B) + C*alpha[:, i])
    

    for i in range(N-1, 0, -1):
        U_new[:, i-1] = alpha[:, i] * U_new[:, i] + beta[:, i]

    # Step 2 calculations
    D = np.zeros_like(U_old)
    alpha = np.zeros_like(U_old)
    beta = np.zeros_like(U_old)

    D[1:-1, 1:-1] = U_new[1:-1, 1:-1]/dt - a * (U_old[2:, 1:-1] - 2*U_old[1:-1, 1:-1] + U_old[:-2, 1:-1])/(2*dy**2)
    beta[1, :] = U_old[0, :]

    for i in range(1, N-1):
        alpha[i+1, :] = -np.full(N, A) / (np.full(N, B) + C*alpha[i, :])
        beta[i+1, :] = (D[i, :] - C*beta[i, :]) / (np.full(N, B) + C*alpha[i, :])
    for i in range(N-1, 0, -1):
        U_old[i-1, :] = alpha[i, :] * U_old[i, :] + beta[i, :]

    # left (x = 0 by y)
    U_new[:M1, 0] = 0 
    U_new[M1:M2, 0] = 1
    U_new[M2:, 0] = 0
    # right
    U_new[:M1, M-1] = 1  
    U_new[M1:, M-1] = 0
    # top
    U_new[N-1, :] = 0  
    # bottom (y = 0 by x)
    U_new[0, M-1] = 0 

    max_diff = np.max(np.abs(U_new - U_old))
    iter += 1

plt.title("Fractional step method")
plt.xlabel("x")
plt.ylabel("y")
plt.contourf(x, y, U_new)
plt.tight_layout
plt.show()