import numpy as np
from numpy import cos, sin, pi
from numpy.linalg import inv
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def graph(T, n):
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    z = np.linspace(0, 1, n)
    x, y, z = np.meshgrid(x, y, z, indexing='ij')

    fig = go.Figure(data=[go.Volume(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=T.flatten(),
        isomin=T.min(),
        isomax=T.max(),
        opacity=0.1,
        surface_count=21,
        colorscale='Viridis')])
    fig.update_layout(scene=dict(
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')))

    fig.show()

start_x, end_x = (0, 1)
start_y, end_y = (0, 1)
start_z, end_z = (0, 1)

N = 21
dx = (end_x - start_x) / (N - 1)
dy = (end_y - start_y) / (N - 1)
dz = (end_z - start_z) / (N - 1)

N1 = int(0.4 * N)
N2 = int(0.7 * N)
M1 = int(0.3 * N)
M2 = int(0.6 * N)
P1 = int(0.7 * N)

x = np.arange(start=0, stop=N) * dx
y = np.arange(start=0, stop=N) * dy
z = np.arange(start=0, stop=N) * dz

U_old = np.zeros((N, N, N))
U_new = np.zeros((N, N, N))

A = np.zeros((N-2, N-2))
B = np.zeros((N-2, N-2))
C = np.zeros((N-2, N-2))
D = np.zeros((N-2, N))

a = np.zeros((N, N))

alpha = np.zeros((N, N-2, N-2))
beta = np.zeros((N-2, N))

# initial condition
U_new[0:N, 0:N, 0:N] = 0

# boundary conditions
U_new[int(0.3 * N):int(0.7 * N), -1, int(0.7 * N):N] = 1
U_new[N-1, M1:M2, 0:int(0.3 * N)] = 1

for j in range(0, N):
    for n in range(1, 100):
        A[:, :] = 0
        B[:, :] = 0
        C[:, :] = 0
        D[:, :] = 0
        beta[:, :] = 0
        alpha[:] = 0
        a[:, :] = 0

        np.fill_diagonal(A, 1 / dx**2)
        np.fill_diagonal(C, 1 / dx**2)
        np.fill_diagonal(B[0:, 1:], 1 / dz**2)
        np.fill_diagonal(B[1:, 0:], 1 / dz**2)
        np.fill_diagonal(B, 2*cos(pi*n/N)/dy**2 - 2/dx**2 - 2/dy**2 - 2/dz**2)

        alpha[1] = 0
        beta[:, 1] = U_new[0, j, 1:N-1]
        for i in range(1, N-1):
            denominator = inv(B + C @ alpha[i])
            alpha[i+1] = -denominator @ A
            beta[:, i+1] = denominator @ (D[:, i] - C @ beta[:, i])

        a[1:N-1, N-1] = U_new[1:N-1, j, N-1]
        for i in range(N-2, -1, -1):
            a[1:N-1, i] = alpha[i+1] @ a[1:N-1, i+1] \
                        + beta[:, i+1]

        U_new[1:N-1, j, 1:N-1] += a[1:N-1, 1:N-1] * sin(pi*j*n/N)

#print(np.max(U_new))
graph(U_new, N)
plt.show()