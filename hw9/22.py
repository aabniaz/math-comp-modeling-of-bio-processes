import numpy as np
import plotly.graph_objects as go

eps = 1e-6
iter_count = 1
h = 0.0447
dt = 0.001

x = np.arange(0, h+1, h)
y = np.arange(0, h+1, h)
z = np.arange(0, h+1, h)
m, n, k = len(x), len(y), len(z)
m1 = int(0.4 * m)
m2 = int(0.6 * m)

A = -1 / (h ** 2)
B = 1 / dt + 2 * 1 / h ** 2
C = -1 / (h ** 2)

u = np.zeros((m, n, k))
u[-1, :m1, m2:] = 1
u[0, m2:, :m1] = 1

u_1 = u.copy()
u_2 = u.copy()

u_old = np.zeros_like(u)  

# Alternating Direction Method
while True:
    # 1
    D = u / dt + (np.roll(u, 1, axis=0) - 2 * u + np.roll(u, -1, axis=0)) / h ** 2 + \
        (np.roll(u, 1, axis=1) - 2 * u + np.roll(u, -1, axis=1)) / h ** 2 + \
        (np.roll(u, 1, axis=2) - 2 * u + np.roll(u, -1, axis=2)) / h ** 2

    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, m - 1):
        alpha[i + 1, :, :] = -np.full(m, A) / (np.full(m, B) + C * alpha[i, :, :])
        beta[i + 1, :, :] = (D[i, :, :] - C * beta[i, :, :]) / (np.full(m, B) + C * alpha[i, :, :])

    for i in range(n - 1, 0, -1):
        u_1[i - 1, :, :] = alpha[i, :, :] * u_1[i, :, :] + beta[i, :, :]

    u_1[-1, :m1, m2:] = 1
    u_1[0, m2:, :m1] = 1

    # 2
    D = u_1 / dt + (np.roll(u_1, 1, axis=0) - 2 * u_1 + np.roll(u_1, -1, axis=0)) / h ** 2 + \
        (np.roll(u_1, 1, axis=1) - 2 * u_1 + np.roll(u_1, -1, axis=1)) / h ** 2 + \
        (np.roll(u_1, 1, axis=2) - 2 * u_1 + np.roll(u_1, -1, axis=2)) / h ** 2

    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, n - 1):
        alpha[:, i + 1, :] = -np.full(n, A) / (np.full(n, B) + C * alpha[:, i, :])
        beta[:, i + 1, :] = (D[:, i, :] - C * beta[:, i, :]) / (np.full(n, B) + C * alpha[:, i, :])

    for i in range(n - 1, 0, -1):
        u_2[:, i - 1, :] = alpha[:, i, :] * u_2[:, i, :] + beta[:, i, :]

    u_2[-1, :m1, m2:] = 1
    u_2[0, m2:, :m1] = 1

    # 3
    D = u_2 / dt + (np.roll(u_2, 1, axis=0) - 2 * u_2 + np.roll(u_2, -1, axis=0)) / h ** 2 + \
        (np.roll(u_2, 1, axis=1) - 2 * u_2 + np.roll(u_2, -1, axis=1)) / h ** 2 + \
        (np.roll(u_2, 1, axis=2) - 2 * u_2 + np.roll(u_2, -1, axis=2)) / h ** 2

    alpha = np.zeros_like(u)
    beta = np.zeros_like(u)
    for i in range(1, k - 1):
        alpha[:, :, i + 1] = -np.full(k, A) / (np.full(k, B) + C * alpha[:, :, i])
        beta[:, :, i + 1] = (D[:, :, i] - C * beta[:, :, i]) / (np.full(k, B) + C * alpha[:, :, i])

    for i in range(k - 1, 0, -1):
        u[:, :, i - 1] = alpha[:, :, i] * u[:, :, i] + beta[:, :, i]

    u[-1, :m1, m2:] = 1
    u[0, m2:, :m1] = 1

    max_diff = np.max(np.abs(u - u_old))
    iter_count += 1
    if max_diff < eps:
        break

    u_old = u.copy()  

print("iterations:", iter_count)

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
