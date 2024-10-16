# import numpy as np
# import matplotlib.pyplot as plt

# a = 1
# eps = 1e-6
# iter = 0
# max = 1
# stop_iter = 10000
# dx = 0.01
# dy = 0.01
# dt = 0.0001

# x = np.arange(0, dx+1, dx)
# y = np.arange(0, dy+1, dy)

# M = len(x)
# N = len(y)

# M1 = int(0.3*M)
# m2 = int(0.7*M)

# A = -a/(2*dx**2)
# B = 1/dt + a/dx**2
# C = -a/(2*dx**2)

# u = np.zeros((M, N))

# def bound_cond(u):
#     # top
#     u[-1, :] = 0

#     # bottom
#     u[0, -1] = 0 #y = 0 by x 

#     # left
#     u[:M1, 0] = 0 #x = 0 by y 
#     u[M1:m2, 0] = 1
#     u[m2:, 0] = 0

#     # right
#     u[:M1, -1] = 1
#     u[M1:, -1] = 0
    
# bound_cond(u)

# u_mid = u.copy()
# while max > eps and iter < stop_iter:
#     # step 1
#     u_old = u.copy()
#     D = np.zeros_like(u)
#     alpha = np.zeros_like(u)
#     beta = np.zeros_like(u)
    
#     D[1:-1, 1:-1] = u[1:-1, 1:-1]/dt + a * (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2])/2*dx**2 + a*(u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/dy**2
    
#     beta[:, 1] = u[:, 0]
    
#     for i in range(1, M-1):
#         alpha[:, i+1] = -np.full(M, A) / (np.full(M, B) + C*alpha[:, i])
#         beta[:, i+1] = (D[:, i] - C*beta[:, i]) / (np.full(M, B) + C*alpha[:, i])
        
#     for i in range(N-1, 0, -1):
#         u_mid[:, i-1] = alpha[:, i] * u_mid[:, i] + beta[:, i]
        
#     bound_cond(u_mid)
#     bound_cond(u)
    
#     # step 2
#     D = np.zeros_like(u)
#     alpha = np.zeros_like(u)
#     beta = np.zeros_like(u)
    
#     D[1:-1, 1:-1] = u_mid[1:-1, 1:-1]/dt - a * (u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/2*dy**2
    
#     beta[1, :] = u[0, :]
    
#     for i in range(1, N-1):
#         alpha[i+1, :] = -np.full(N, A) / (np.full(N, B) + C*alpha[i, :])
#         beta[i+1, :] = (D[i, :] - C*beta[i, :]) / (np.full(N, B) + C*alpha[i, :])
        
#     for i in range(N-1, 0, -1):
#         u[i-1, :] = alpha[i, :] * u[i, :] + beta[i, :]
        
#     bound_cond(u)
    
#     max_d = np.max(np.abs(u - u_old))
    
#     iter += 1
    
#     if max_d < eps:
#         break
    
# print(iter)
# plt.title("Fractional step method")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.contourf(x, y, u)
# plt.tight_layout
# plt.show()


