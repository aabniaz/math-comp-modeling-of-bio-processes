#Gauss-Seidel method
import numpy as np
import matplotlib.pyplot as plt

iter = 1
eps = 1e-6
deltax = 0.01
deltay = 0.01

x = np.arange(0, 1+deltax, deltax)
y = np.arange(0, 1+deltay, deltay)

M = len(x)
N = len(y)
M1 = int(0.3/deltay)
M2 = int(0.7/deltay)

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

arr = [U]

while True:
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

    for i in range(1, N-1):
        for j in range(1, M-1):
            U_new[i,j] = (1/deltax**2 * (arr[-1][i, j+1] + U_new[i, j-1]) + 1/deltay**2 * (arr[-1][i+1, j] + U_new[i-1, j])) / (2/deltax**2 + 2/deltay**2)
    arr.append(U_new)
    
    max_diff = np.max(np.abs(U - U_new))
    print(max_diff)
    
    U[:, :] = U_new[:, :]
    iter += 1
             
    if eps > max_diff:
        break
        
print(iter)
plt.contourf(x, y, U)
plt.show()