import numpy as np
import matplotlib.pyplot as plt

iter = 0
dx = 0.001 #0.01
# eps = 0.1

def analytical(x):
    return -x**4 / 12 - x**3 / 3 - 3 * x**2 / 2 + 11 * x / 12 + 1
x = np.arange(0, dx+1, dx)
n = len(x)

P_a = analytical(x)
print('Analytical:\n', P_a)

def func(x):
    return -x**2 - 2*x - 3

P = np.zeros(n)
# b. c
P[0] = 1
P[-1] = 0

A = -1/(12*dx**2)
B = 16/(12*dx**2)
C = -30/(12*dx**2)
D = 16/(12*dx**2)
E = -1/(12*dx**2)
H = func(x)

alpha = np.zeros(n)
beta = np.zeros(n)
gamma = np.zeros(n)
alpha[0] = 0  # from b.c
beta[0] = 0
gamma[0] = 1

alpha[1] = -(B + D*beta[0]) / (C + D*alpha[0])
beta[1] = -A / (C + D*alpha[0])
gamma[1] = (H[0] - D*gamma[0]) / (C + D*alpha[0])

# while True:
for i in range(1, n-1):
    alpha[i+1] = - ((B + D*beta[i] + E*alpha[i-1]*beta[i])/ (C + D*alpha[i] + E*alpha[i-1]*alpha[i] + E*beta[i-1]))
    beta[i+1] = - (A / (C + D*alpha[i] + E*alpha[i-1]*alpha[i] + E*beta[i-1]))
    gamma[i+1] = (H[i] - D*gamma[i] - E*alpha[i-1]*gamma[i] - E*gamma[i-1]) / (C + D*alpha[i] + E*alpha[i-1]*alpha[i] + E*beta[i-1])


P[-2] = alpha[-1]*P[-1] + beta[-1]*P[0] + gamma[-1]

for i in range(n-3, 0, -1):
    P[i] = alpha[i+1]*P[i+1] + beta[i+1]*P[i+2] + gamma[i+1]

#     max_diff = np.max(np.abs(P[i] - P[i+1]))

#     if max_diff < eps:
#         break

# iter += 1    
# print("Number of iterations: ", iter)

print('\nFive diagonal numerical:\n', P)
# print(H)
print("Max Error ", np.max(np.abs(P - P_a)))
    
plt.plot(x, P, label = 'numerical')
plt.plot(x, P_a, label = 'analytical')
plt.legend()
plt.show()