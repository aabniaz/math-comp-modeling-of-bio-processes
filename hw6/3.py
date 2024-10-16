#Relaxation method
import numpy as np
import matplotlib.pyplot as plt

omega = 1.2
h = 0.01
iter = 1
eps = 1e-6

x = np.arange(0, h+1, h)
y = np.arange(0, h+1, h)

M = len(x)
N = len(y)
M1 = int(0.3/h)
M2 = int(0.6/h)
U = np.zeros((N, M))

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

while(True):
    iter += 1
    U_new = np.zeros((N, M))
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
            U_new[i,j] = omega/4*((arr[-1][i+1,j] + U_new[i-1,j] + arr[-1][i,j+1] + U_new[i,j-1]) - 4*(1 - 1/omega)*arr[-1][i,j])
            
    arr.append(U_new)
    
    max_diff = np.max(np.abs(arr[-1] - arr[-2]))
    print(max_diff)
    
    U[:, :] = U_new[:, :]
    
    if max_diff < eps:
        break

print(iter)
plt.contourf(x, y, U_new)
plt.show()