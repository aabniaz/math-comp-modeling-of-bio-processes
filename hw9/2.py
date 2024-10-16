import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

a = 1
eps = 1e-6
iter = 1
h = 0.01
dt = 0.00001
x = np.arange(0, h+1, h)
y = np.arange(0, h+1, h)
z = np.arange(0, h+1, h)
m, n, k = len(x), len(y), len(z)
m1 = int(0.4*m)
m2 = int(0.6*m)
A = -a/(h**2)
B = 1/dt + 2*a/h**2
C = -a/(h**2)
u = np.zeros((m, n, k))
def bound_cond(u):
# x y z
    u[-1, :m1, m2:] = 1
    u[0, m2:, :m1] = 1
bound_cond(u)
u_1 = u.copy()
u_2 = u.copy()
while True:
# step 1
    u_old = u.copy()
    D = np.zeros_like(u)
    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    D[1:-1, 1:-1, 1:-1] = u[1:-1, 1:-1, 1:-1]/dt + (u[1:-1, 2:, 1:-1] - 2*u[1:-1, 1:-1, 
    1:-1] + u[1:-1, :-2, 1:-1])/h**2 + (u[1:-1, 1:-1, 2:] - 2*u[1:-1, 1:-1, 1:-1] + u[1:-1, 
    1:-1, :-2])/h**2
    beta[1, :, :] = u[0, :, :]
    for i in range(1, m-1):
        alpha[i+1, :, :] = -np.full(m, A) / (np.full(m, B) + C*alpha[i, :, :])
        beta[i+1, :, :] = (D[i, :, :] - C*beta[i, :, :]) / (np.full(m, B) + C*alpha[i, :, :])
    for i in range(n-1, 0, -1):
        u_1[i-1, :, :] = alpha[i, :, :] * u_1[i, :, :] + beta[i, :, :]
    bound_cond(u_1)
    bound_cond(u)
# step 2
    D = np.zeros_like(u)
    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    D[1:-1, 1:-1, 1:-1] = u_1[1:-1, 1:-1, 1:-1]/dt + (u_1[2:, 1:-1, 1:-1] - 2*u_1[1:-1, 
    1:-1, 1:-1] + u_1[:-2, 1:-1, 1:-1])/h**2 + (u_1[1:-1, 1:-1, 2:] - 2*u_1[1:-1, 1:-1, 1:-1] 
    + u_1[1:-1, 1:-1, :-2])/h**2
    beta[:, 1, :] = u[:, 0, :]
    for i in range(1, n-1):
        alpha[:, i+1, :] = -np.full(n, A) / (np.full(n, B) + C*alpha[:, i, :])
        beta[:, i+1, :] = (D[:, i, :] - C*beta[:, i, :]) / (np.full(n, B) + C*alpha[:, i, :])
    for i in range(n-1, 0, -1):
        u_2[:, i-1, :] = alpha[:, i, :] * u_2[:, i, :] + beta[:, i, :]
    bound_cond(u_2)
    bound_cond(u)
# step 3
    D = np.zeros_like(u)
    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    D[1:-1, 1:-1, 1:-1] = u_2[1:-1, 1:-1, 1:-1]/dt + (u_2[2:, 1:-1, 1:-1] - 2*u_2[1:-1, 
    1:-1, 1:-1] + u_2[:-2, 1:-1, 1:-1])/h**2 + (u_2[1:-1, 2:, 1:-1] - 2*u_2[1:-1, 1:-1, 1:-1] 
    + u_2[1:-1, :-2, 1:-1])/h**2
    beta[:, :, 1] = u[:, :, 0]
    for i in range(1, k-1):
        alpha[:, :, i+1] = -np.full(k, A) / (np.full(k, B) + C*alpha[:, :, i]) 
        beta[:, :, i+1] = (D[:, :, i] - C*beta[:, :, i]) / (np.full(k, B) + C*alpha[:, :, i])
    for i in range(k-1, 0, -1):
        u[:, :, i-1] = alpha[:, :, i] * u[:, :, i] + beta[:, :, i]
    bound_cond(u)
    max_d = np.max(np.abs(u - u_old))
    iter += 1
    if max_d < eps:
        break

print(iter)
x, y, z = np.meshgrid(x, y, z, indexing='ij')
fig = go.Figure(data=[go.Volume(
x=x.flatten(),
y=y.flatten(),
z=z.flatten(),
value=u.flatten(),
isomin=0,
isomax=1,
opacity=0.1, 
surface_count=21, 
colorscale='Viridis')])
fig.update_layout(scene=dict(
xaxis=dict(title='X'),
yaxis=dict(title='Y'),
zaxis=dict(title='Z')))
fig.show()