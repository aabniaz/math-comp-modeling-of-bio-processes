import numpy as np 
import matplotlib.pyplot as plt 
 
eps = 1e-6 
iter = 0 
Re = 100 
stop_iter = 10000 
 
dt = 0.01 
dx = 0.01 
dy = 0.01 
 
t = np.arange(0, dt+1, dt) 
x = np.arange(0, dx+1, dx) 
y = np.arange(0, dy+1, dy) 
 
m, n, p = len(x), len(y), len(t) 
 
m1 = int(0.3*m) 
m2 = int(0.6*m) 
 
u = np.zeros((m, n)) 
v = np.zeros((m, n)) 
 
def bound_u(u): 
    # inlet  
    #   y     x 
    u[-1, m2:] = 0 
     
    # outlet 
    u[m1:m2, 0] = u[m1:m2, 1] 
    u[0, m2:] = u[1, m2:] 
     
    # walls 
    u[0, :m2] = 0  
    u[:, -1] = 0 
    u[:m1, 0] = 0 
    u[m2:, 0] = 0 
    u[-1, :m2] = 0 
     
def bound_v(v): 
    # inlet  
    #   y     x 
    v[-1, m2:] = -1 
     
    # outlet 
    v[m1:m2, 0] = v[m1:m2, 1] 
    v[0, m2:] = v[1, m2:] 
     
    # walls 
    v[0, :m2] = 0 
    v[m2:, 0] = 0 
    v[:, -1] = 0
    v[:m1, 0] = 0 
    v[-1, :m2] = 0 
     
bound_u(u) 
bound_v(v) 
 
u_05 = u.copy() 
v_05 = v.copy() 
 
# bugers eq solving: 
while True: 
    u_old = u.copy() 
    v_old = v.copy() 
     
    # u 
    # step 1:  
    A = np.zeros_like(u) 
    B = np.zeros_like(u) 
    C = np.zeros_like(u) 
    D = np.zeros_like(u) 
 
    A[1:-1, 1:-1] = -0.5*(1/(Re*dx**2) - u[1:-1, 1:-1]/dx) 
    B[1:-1, 1:-1] = 1/dt - 0.5*(u[1:-1, 1:-1]/dx - 2/(Re*dx**2)) 
    C[1:-1, 1:-1] = -0.5/(Re*dx**2) 
    D[1:-1, 1:-1] = u[1:-1, 1:-1]/dt + 0.5*(1/Re * (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2])/dx**2 - u[1:-1, 1:-1]*(u[1:-1, 2:] - u[1:-1, 1:-1])/dx) + 1/Re*(u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/dy**2 - v[1:-1, 1:-1]*(u[2:, 1:-1] - u[1:-1, 1:-1])/dy 
 
    alpha = np.zeros_like(u) 
    beta = np.zeros_like(u) 
 
    alpha[:, 1] = 0 
    beta[:, 1] = v[:, 0] 
 
    for i in range(1, m-1): 
        alpha[1:-1, i+1] = -A[1:-1, i] / (B[1:-1, i] + alpha[1:-1, i]*C[1:-1, i]) 
        beta[1:-1, i+1] = (D[1:-1, i] - C[1:-1, i]*beta[1:-1, i]) / (B[1:-1, i] + alpha[1:-1, i]*C[1:-1, i]) 
 
    for i in range(m-1, 0, -1): 
        u_05[:, i-1] = alpha[:, i]*u_05[:, i] + beta[:, i] 
 
    bound_u(u_05) 
    bound_u(u) 
     
    # u 
    # step 2 
    A = np.zeros_like(u) 
    B = np.zeros_like(u) 
    C = np.zeros_like(u) 
    D = np.zeros_like(u) 
 
    A[1:-1, 1:-1] = -0.5 * ((1/Re) * (1/dy**2) - v[1:-1, 1:-1]/dy) 
    B[1:-1, 1:-1] = 1/dt - 0.5*(v[1:-1, 1:-1]/dy - (2/Re) * (1/dy**2)) 
    C[1:-1, 1:-1] = -0.5/(Re*dy**2) 
    D[1:-1, 1:-1] = u_05[1:-1, 1:-1]/dt + 0.5*(v[1:-1, 1:-1]*(u[2:, 1:-1] - u[1:-1, 1:-1])/dy - 1/Re*(u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1])/dy**2) 
 
    alpha = np.zeros_like(u) 
    beta = np.zeros_like(u) 
 
    alpha[1, :] = 0 
    beta[1, :] = v[1, :] 
 
    for i in range(1, n-1): 
        alpha[i+1, 1:-1] = -A[i, 1:-1] / (B[i, 1:-1] + alpha[i, 1:-1]*C[i, 1:-1]) 
        beta[i+1, 1:-1] = (D[i, 1:-1] - C[i, 1:-1]*beta[i, 1:-1]) / (B[i, 1:-1] + alpha[i, 1:-1]*C[i, 1:-1]) 
 
    for i in range(n-1, 0, -1): 
        u[i-1, :] = alpha[i, :]*u[i, :] + beta[i, :] 
 
    bound_u(u) 
 
    max_u = np.max(np.abs(u - u_old)) 
     
     
    # v    
    # step 1 
    A = np.zeros_like(v) 
    B = np.zeros_like(v) 
    C = np.zeros_like(v) 
    D = np.zeros_like(v) 
 
    A[1:-1, 1:-1] = -0.5*(1/(Re*dx**2) - u[1:-1, 1:-1]/dx) 
    B[1:-1, 1:-1] = 1/dt - 0.5*(u[1:-1, 1:-1]/dx - 2/(Re*dx**2)) 
    C[1:-1, 1:-1] = -0.5/(Re*dx**2) 
    D[1:-1, 1:-1] = v[1:-1, 1:-1]/dt + 0.5*(1/Re * (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, :-2])/dx**2 - u[1:-1, 1:-1]*(v[1:-1, 2:] - v[1:-1, 1:-1])/dx) + 1/Re*(v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1])/dy**2 - v[1:-1, 1:-1]*(v[2:, 1:-1] - v[1:-1, 1:-1])/dy 
 
    alpha = np.zeros_like(v) 
    beta = np.zeros_like(v) 
 
    alpha[:, 1] = 0 
    beta[:, 1] = v[:, 0] 
     
    for i in range(1, m-1): 
        alpha[1:-1, i+1] = -A[1:-1, i] / (B[1:-1, i] + alpha[1:-1, i]*C[1:-1, i]) 
        beta[1:-1, i+1] = (D[1:-1, i] - C[1:-1, i]*beta[1:-1, i]) / (B[1:-1, i] + alpha[1:-1, i]*C[1:-1, i]) 
 
    for i in range(m-1, 0, -1): 
        v_05[:, i-1] = alpha[:, i]*v_05[:, i] + beta[:, i] 
 
    bound_v(v_05) 
    bound_v(v) 
 
    #v 
    # step 2 
    A = np.zeros_like(v) 
    B = np.zeros_like(v) 
    C = np.zeros_like(v) 
    D = np.zeros_like(v) 
 
    A[1:-1, 1:-1] = -0.5 * ((1/Re) * (1/dy**2) - v[1:-1, 1:-1]/dy) 
    B[1:-1, 1:-1] = 1/dt - 0.5*(v[1:-1, 1:-1]/dy - (2/Re) * (1/dy**2)) 
    C[1:-1, 1:-1] = -0.5/(Re*dy**2) 
    D[1:-1, 1:-1] = v_05[1:-1, 1:-1]/dt + 0.5*(v[1:-1, 1:-1]*(v[2:, 1:-1] - v[1:-1, 1:-1])/dy - 1/Re*(v[2:, 1:-1] - 2*v[1:-1, 1:-1] + v[:-2, 1:-1])/dy**2) 
 
    alpha = np.zeros_like(v) 
    beta = np.zeros_like(v) 
     
    alpha[1, :] = 0 
    beta[1, :] = v[1, :] 
 
    for i in range(1, n-1): 
        alpha[i+1, 1:-1] = -A[i, 1:-1] / (B[i, 1:-1] + alpha[i, 1:-1]*C[i, 1:-1]) 
        beta[i+1, 1:-1] = (D[i, 1:-1] - C[i, 1:-1]*beta[i, 1:-1]) / (B[i, 1:-1] + alpha[i, 1:-1]*C[i, 1:-1]) 
 
    for i in range(n-1, 0, -1): 
        v[i-1, :] = alpha[i, :]*v[i, :] + beta[i, :] 
 
    bound_v(v) 
 
    max_v = np.max(np.abs(v - v_old)) 
 
    iter += 1 
     
    if (max_u < eps and max_v < eps) or iter > stop_iter: 
        break 
 
print(iter) 
 
X, Y = np.meshgrid(x, y) 
plt.contourf(X, Y, np.sqrt(u**2 + v**2), levels = 10)  
plt.title('u + v') 
plt.show()