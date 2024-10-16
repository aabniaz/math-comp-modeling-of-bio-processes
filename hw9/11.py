import numpy as np
import plotly.graph_objects as go

a = 1
eps = 1e-6
max_iter = 1
h = 0.01
dt = 0.0001

x = np.arange(0, h+1, h)
y = np.arange(0, h+1, h)
z = np.arange(0, h+1, h)

m, n, k = len(x), len(y), len(z)
m1 = int(0.3*m)
m2 = int(0.7*m)

A = -a/(2*h**2)
B = 1/dt + a/h**2
C = -a/(2*h**2) 

u = np.zeros((m, n, k))

u[m1:m2, -1, m2:] = 1
u[-1, m1:m2, :m1] = 1

u_old = u.copy()

# Fractional Step Method
while True:
    # 1
    D = u_old/dt + a*0.5 * ((np.roll(u_old, 1, axis=0) - 2*u_old + np.roll(u_old, -1, axis=0))/h**2 +
                         (np.roll(u_old, 1, axis=1) - 2*u_old + np.roll(u_old, -1, axis=1))/h**2 +
                         (np.roll(u_old, 1, axis=2) - 2*u_old + np.roll(u_old, -1, axis=2))/h**2)

    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, m-1):
        alpha[i+1, :, :] = -A / (B + C*alpha[i, :, :])
        beta[i+1, :, :] = (D[i, :, :] - C*beta[i, :, :]) / (B + C*alpha[i, :, :])
        
    u_1 = np.zeros_like(u)
    for i in range(n-1, 0, -1):
        u_1[i-1, :, :] = alpha[i, :, :] * u_1[i, :, :] + beta[i, :, :]

    u_1[m1:m2, -1, m2:] = 1
    u_1[-1, m1:m2, :m1] = 1

    # 2
    D = u_1/dt - a*0.5 * ((np.roll(u, -1, axis=1) - 2*u + np.roll(u, 1, axis=1))/h**2)
    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, n-1):
        alpha[:, i+1, :] = -A / (B + C*alpha[:, i, :])
        beta[:, i+1, :] = (D[:, i, :] - C*beta[:, i, :]) / (B + C*alpha[:, i, :])
        
    u_2 = np.zeros_like(u)
    for i in range(n-1, 0, -1):
        u_2[:, i-1, :] = alpha[:, i, :] * u_2[:, i, :] + beta[:, i, :]

    u_2[m1:m2, -1, m2:] = 1
    u_2[-1, m1:m2, :m1] = 1

    # 3
    D = u_2/dt - a*0.5 * ((np.roll(u, -1, axis=2) - 2*u + np.roll(u, 1, axis=2))/h**2)
    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, k-1):
        alpha[:, :, i+1] = -A / (B + C*alpha[:, :, i])
        beta[:, :, i+1] = (D[:, :, i] - C*beta[:, :, i]) / (B + C*alpha[:, :, i])
        
    u_new = np.zeros_like(u)
    for i in range(k-1, 0, -1):
        u_new[:, :, i-1] = alpha[:, :, i] * u_new[:, :, i] + beta[:, :, i]

    u_new[m1:m2, -1, m2:] = 1
    u_new[-1, m1:m2, :m1] = 1

    max_diff = np.max(np.abs(u_new - u_old))
    max_iter += 1
    if max_diff < eps:
        break
    u_old = u_new.copy()

print("Iterations:", max_iter)

x, y, z = np.meshgrid(x, y, z, indexing='ij')
fig = go.Figure(data=[go.Volume(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value=u_new.flatten(),
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