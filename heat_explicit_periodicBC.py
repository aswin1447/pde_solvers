import numpy as np
import matplotlib.pyplot as plt

L = 1.
t = .1
omega = 2.

steps = 1e3
dt = t / steps
N = 100.  # Try N = 150.
dx = L / N

D = .1
A = -1.
B = .001

# b<=1 !!!
b = 2. * D * dt / (dx ** 2)

print(dt, dx, b)

X = np.linspace(0, L, N + 1)
T = np.linspace(0, t, steps + 1)

# initial temperature distribution
u = np.sin(2 * np.pi * omega * X)
u = (1 / np.sqrt(B)) * np.ones(len(X))
u = np.zeros(len(X))
u[10:20] = 1.

u = np.array(list(map(lambda x: u, T)))

# for k in range(0, int(steps)):
#     for l in range(0, int(N + 1)):
#         if l - 1 == -1:
#             print(k, l, 'Here!')
#             u[k + 1][l] = u[k][l] + 0.5 * b * (u[k][-1] + u[k][l + 1] - 2. * u[k][l])
#         elif l + 1 == N + 1:
#             print(k, l, 'Here!')
#             u[k + 1][l] = u[k][l] + 0.5 * b * (u[k][l - 1] + u[k][0] - 2. * u[k][l])
#         else:
#             print(k, l, 'Here!')
#             u[k + 1][l] = u[k][l] + 0.5 * b * (u[k][l - 1] + u[k][l + 1] - 2. * u[k][l])
#
#         u[k + 1][l] += A * u[k][l] + B * (u[k][l] ** 3)

for k in range(0, int(steps)):
    for l in range(1, int(N)):
        u[k + 1][l] = u[k][l] + 0.5 * b * (u[k][l - 1] + u[k][l + 1] - 2. * u[k][l])
        # u[k + 1][l] += A * u[k][l] + B * (u[k][l] ** 3)
        # print(A * u[k][l] + B * (u[k][l] ** 3))

filtered = int(steps / 100)

plt.plot(X, u[0])
# plt.plot(X, u[5])
# plt.plot(X, u[10])
# plt.plot(X, u[32])
plt.plot(X, u[100])
plt.plot(X, u[200])
plt.plot(X, u[900])
plt.plot(X, u[int(steps)])
# plt.plot(X, np.exp(- D * np.pi ** 2 * t) * np.sin(np.pi * X), color='k', alpha=.5)

# x, t = np.meshgrid(X, T[0::filtered])
# plt.pcolormesh(x, t, u[0::filtered])
# plt.colorbar()

plt.ticklabel_format(useOffset=False)
plt.xlim(0, 1)
plt.tight_layout()
plt.savefig('plots/explicit_heat_periodicBC.pdf')
