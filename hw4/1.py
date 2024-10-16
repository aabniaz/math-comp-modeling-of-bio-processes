import numpy as np
import matplotlib.pyplot as plt

deltax = 0.001 
x = np.arange(0, 1 + deltax, deltax)  
N = len(x)

A = -1 / (12 * deltax**2)
B = 16 / (12 * deltax**2)
C = -30 / (12 * deltax**2)
D = 16 / (12 * deltax**2)
E = -1 / (12 * deltax**2)
H = -x**2 - 2*x - 3

def five_diagonal_matrix():
    alpha = np.zeros(N)
    beta = np.zeros(N)
    gamma = np.zeros(N)
    alpha[0] = 0 
    beta[0] = 0
    gamma[0] = 1

    alpha[1] = -(B + D * beta[0]) / (C + D * alpha[0])
    beta[1] = -A / (C + D * alpha[0])
    gamma[1] = (H[0] - D * gamma[0]) / (C + D * alpha[0])

    for i in range(1, N-1):
        denominator = C + D * alpha[i] + E * alpha[i-1] * alpha[i] + E * beta[i-1]
        alpha[i+1] = -(B + D * beta[i] + E * alpha[i-1] * beta[i]) / denominator
        beta[i+1] = -A / denominator
        gamma[i+1] = (H[i] - D * gamma[i] - E * alpha[i-1] * gamma[i] - E * gamma[i-1]) / denominator

    P = np.zeros(N)
    P[-2] = (H[-2] - D * gamma[-2] - E * alpha[-3] * gamma[-2] - E * gamma[-3]) / (C + D * alpha[-2] + E * alpha[-3] * alpha[-2] + E * beta[-3])

    for i in range(N-3, 0, -1):
        P[i] = alpha[i+1] * P[i+1] + beta[i+1] * P[i+2] + gamma[i+1]

    P[0] = 1
    P[-1] = 0
    return x, P

def analytical_sol(x):
    return -x**4 / 12 - x**3 / 3 - 3 * x**2 / 2 + 11 * x / 12 + 1

x, numerical_solution = five_diagonal_matrix()
analytical_solution = analytical_sol(x)
max_error = np.max(np.abs(numerical_solution - analytical_solution))
print(f"maximum error: {max_error}")

plt.plot(x, numerical_solution, label='Numerical Solution')
plt.plot(x, analytical_solution, label='Analytical Solution')
plt.xlabel('x')
plt.ylabel('P(x)')
plt.legend()
plt.show()