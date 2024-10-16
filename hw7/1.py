import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

dx = 0.01
dy = 0.01

x = int(1/dx)
y = int(1/dy)
n = x

P = np.zeros((n, n))
alpha = np.zeros((n, n, n))
betta = np.zeros((n, n))

A = np.zeros((n, n))
B = np.zeros((n, n))
C = np.zeros((n, n))
D = np.zeros((n, n))
I_min = np.zeros((n, n))

for i in range(0,n):
    A[i][i] = 1/(dx**2)
    C[i][i] = 1/(dx**2)
    B[i][i] = (-2/(dx**2))-(2/(dy**2))
    I_min[i][i]=-1
    if(i<n-1):
        B[i][i+1] = 1/(dy**2)
        B[i+1][i] = 1/(dy**2)

for i in range(33):
    betta[0][i+33] = 1
    P.T[n-1][i] = 1


for i in range(0, n-1):
    alpha[i+1] = np.matmul(I_min, np.matmul(inv(B+np.matmul(C, alpha[i])), A))
    betta[i+1] = np.matmul(inv(B+np.matmul(C, alpha[i])), (D[i]-np.matmul(C,betta[i])))

for i in range(n-1, 0, -1):
    P.T[i-1]=np.matmul(alpha[i], P.T[i])+betta[i]
    
#print(P)
plt.contourf(P)
plt.show()