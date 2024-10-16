# import numpy as np
# import matplotlib.pyplot as plt

# def thomas_algorithm(A, B, C, D):
#     N = len(D)
    
#     alpha = np.zeros(N - 1)
#     beta = np.zeros(N - 1)

#     alpha[0] = -A[0] / B[0]
#     beta[0] = D[0] / B[0]

#     for i in range(1, N-1):
#         denom = B[i] + C[i-1] * alpha[i-1]
#         alpha[i] = -A[i] / denom
#         beta[i] = (D[i] - C[i-1] * beta[i-1]) / denom

#     U = np.zeros(N)
#     U[-1] = (D[-1] - C[-2] * beta[-1]) / (B[-1] + C[-2] * alpha[-1])

#     for i in range(N-2, -1, -1):
#         U[i] = alpha[i] * U[i+1] + beta[i]

#     return U

# N = 100 
# deltax = 1 / (N - 1)
# b = 1  # Heat conduction coefficient
# deltat = 0.01  
# time_steps = 100

# A = -b**2 / deltax**2
# B = 2 * b**2 / deltax**2
# C = -b**2 / deltax**2

# U_0 = np.sin(np.pi * np.linspace(0, 1, N))

# U_current = U_0.copy()
# for _ in range(time_steps):
#     D = U_current / deltat
#     U_next = thomas_algorithm(A, B, C, D)
#     U_current = U_next

# plt.plot(np.linspace(0, 1, N), U_next, label='Numerical Solution')
# plt.xlabel('x')
# plt.ylabel('Temperature (u)')
# plt.legend()
# plt.show()



import numpy as np
import matplotlib.pyplot as plt

def thomas_algorithm(A, B, C, D):
    N = len(D)
    
    # Forward sweep
    alpha = np.zeros(N - 1)
    beta = np.zeros(N - 1)

    alpha[0] = -A[0] / B[0]
    beta[0] = D[0] / B[0]

    for i in range(1, N-1):
        denom = B[i] + C[i-1] * alpha[i-1]
        alpha[i] = -A[i] / denom
        beta[i] = (D[i] - C[i-1] * beta[i-1]) / denom

    # Backward sweep
    U = np.zeros(N)
    U[-1] = (D[-1] - C[-2] * beta[-1]) / (B[-1] + C[-2] * alpha[-1])

    for i in range(N-2, -1, -1):
        U[i] = alpha[i] * U[i+1] + beta[i]

    return U

# Define grid and coefficients
N = 100
deltax = 1 / (N - 1)
deltat = 0.01
x = np.linspace(0, 1, N)
A = -1 / deltax**2 * np.ones(N - 1)
B = 1 / deltat + 2 / deltax**2 * np.ones(N)
C = -1 / deltax**2 * np.ones(N - 1)
D = np.zeros(N)
D[0] = 0  # Initial condition: U(0, t) = 0

# Initialize U
U = np.zeros(N)

# Apply Thomas algorithm for time-stepping
time_steps = 100
for _ in range(time_steps):
    D[1:-1] = U[1:-1] / deltat  # Update D based on previous solution
    U = thomas_algorithm(A, B, C, D)

# Plot the final solution
plt.plot(x, U, label='Numerical Solution')
plt.xlabel('x')
plt.ylabel('U(x, t)')
plt.legend()
plt.show()
