import numpy as np
import matplotlib.pyplot as plt

c = 1
deltat = 0.01
deltax = 0.1
tmax = 1.0
xmax = 1.0

if c * deltat / deltax >= 1:
    raise ValueError(f'cfl condition not satisfied')
print(f'cfl: {c * deltat / deltax}')

timesteps = int(tmax / deltat)
xpoints = int(xmax / deltax)

# init cond
u = np.zeros((timesteps + 1, xpoints + 1))
x = np.linspace(0, xmax, xpoints + 1)
u[0, :] = np.cos(np.pi * x / 2)

# bound cond
u[:, 0] = 1  # u(t, 0) = 1
u[:, -1] = 0  # u(t, 1) = 0

for n in range(timesteps):
    for i in range(1, xpoints + 1):
        u[n + 1, i] = u[n, i] - c * deltat / deltax * (u[n, i] - u[n, i - 1])

for n in range(0, timesteps + 1, int(timesteps / 10)):  
    plt.plot(x, u[n, :], label=f't = {n * deltat:.2f}')

plt.title('numerical solution of transport equation')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.legend()
plt.show()




# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

# c = 1
# deltat = 0.01
# deltax = 0.1
# tmax = 1.0
# xmax = 1.0

# if c * deltat / deltax >= 1:
#     raise ValueError(f'cfl condition not satisfied')
# print(f'cfl: {c * deltat / deltax}')

# timesteps = int(tmax / deltat)
# xpoints = int(xmax / deltax)

# # init cond
# u = np.zeros((timesteps + 1, xpoints + 1))
# x = np.linspace(0, xmax, xpoints + 1)
# u[0, :] = np.cos(np.pi * x / 2)

# # bound cond
# u[:, 0] = 1  # u(t, 0) = 1
# u[:, -1] = 0  # u(t, 1) = 0

# fig, ax = plt.subplots()
# lines, = ax.plot(x, u[0, :])
# ax.set_title('numerical solution of transport equation')
# ax.set_xlabel('x')
# ax.set_ylabel('u(x, t)')
# line_label = ax.text(0.8, 0.9, f't = 0.00', transform=ax.transAxes)

# def update(n):
#     lines.set_ydata(u[n, :])
#     line_label.set_text(f't = {n * deltat:.2f}')
#     for i in range(1, xpoints + 1):
#         u[n + 1, i] = u[n, i] - c * deltat / deltax * (u[n, i] - u[n, i - 1])
#     return lines, line_label

# animation = FuncAnimation(fig, update, frames=timesteps, interval=50, blit=True)
# animation.save('hw5.gif', writer='imagemagick')
# plt.show()