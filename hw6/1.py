#Jacobi method
import numpy as np
import matplotlib.pyplot as plt

iter = 1
eps = 1e-6
deltax = 0.01
deltay = 0.01

x = np.arange(0, deltax + 1, deltax)
y = np.arange(0, deltay + 1, deltay)

M = len(x)
N = len(y)
M1 = int(0.3 / deltay)
M2 = int(0.6 / deltay)

U = np.zeros((N, M))
U_new = np.zeros((N, M))

# left
U[:M1, 0] = 0 #x = 0 by y 
U[M1:M2, 0] = 1
U[M2:, 0] = 0
# right
U[:M1, M-1] = 1
U[M1:, M-1] = 0
# bottom
U[0, M-1] = 0 #y = 0 by x 
# top
U[N-1, :] = 0

while True:
    U_new[1:-1, 1:-1] = (1/deltay**2 * (U[1:-1, 2:] + U[1:-1, 0:-2]) + 1/deltax**2 * (U[2:, 1:-1] + U[0:-2, 1:-1])) / (2/deltax**2 + 2/deltay**2)
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
    max_diff = np.max(np.abs(U - U_new))
    U[:, :] = U_new[:, :]
    iter += 1
    if eps > max_diff:
        break

print(iter)
plt.contourf(x, y, U)
plt.show()