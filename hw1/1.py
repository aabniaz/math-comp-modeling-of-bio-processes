import numpy as np
import matplotlib.pyplot as plt

def thomas_algorithm(A, B, C, D):
    N = len(D)
    alpha = np.zeros(N - 1)
    beta = np.zeros(N - 1)

    alpha[0] = 0
    beta[0] = 1
    
    for i in range(1, N-1):
        denominator = B[i] + C[i-1] * alpha[i-1]
        alpha[i] = -A[i] / denominator
        beta[i] = (D[i] - C[i-1] * beta[i-1]) / denominator

    P = np.zeros(N)
    P[-1] = (D[-1] - C[-2] * beta[-1]) / (B[-1] + C[-2] * alpha[-1])

    for i in range(N-2, 0, -1): 
        P[i] = alpha[i] * P[i+1] + beta[i] 

    P[0]=1
    return P, alpha, beta

N = 400
deltax = 1 / (N - 1)
x = np.linspace(0, 1, N)
A = np.ones(N - 1) / deltax**2
B = -2 * np.ones(N) / deltax**2
C = np.ones(N - 1) / deltax**2
D = -x**2 - 2*x - 3

numerical_solution, alpha, beta = thomas_algorithm(A, B, C, D)

analytical_solution = -x**4 / 12 - x**3 / 3 - 3 * x**2 / 2 + 11 * x / 12 + 1
max_error = np.max(np.abs(numerical_solution - analytical_solution))
appr_error = np.sum(np.abs(numerical_solution - analytical_solution) * deltax)
max_diff = np.max(np.abs(alpha - beta))
num_iter = len(alpha)

print(f"maximum error: {max_error}")
print(f"approximation error: {appr_error}")
print(f"maximum difference: {max_diff}")
print(f"number of iterations: {num_iter}")

plt.plot(x, numerical_solution, label='Numerical Solution')
plt.plot(x, analytical_solution, label='Analytical Solution', linestyle='--')
plt.xlabel('x')
plt.ylabel('P(x)')
plt.legend()
plt.show()