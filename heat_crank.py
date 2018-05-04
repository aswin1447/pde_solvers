import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg

L = 1.
t = .1
omega = 2.

steps = 1e3
dt = t / steps
N = 100.
dx = L / N
D = 0.25

alpha = D * dt / (dx ** 2)

print(dt, dx, 2 * alpha)

X = np.linspace(0, L, N + 1)
T = np.linspace(0, t, steps + 1)

# initial temperature distribution
u = np.sin(np.pi * X)
u = np.cos(2 * np.pi * omega * X) ** 2

# left matrix
mainDiag = 2. * alpha * np.ones(int(N - 1))
mainDiag += 2 * np.ones(int(N - 1))
mainDiag = np.insert(mainDiag, 0, 1)
mainDiag = np.insert(mainDiag, len(mainDiag), 1)
upperDiag = np.insert(-alpha * np.ones(int(N - 1)), 0, 0)
lowerDiag = np.insert(-alpha * np.ones(int(N - 1)), int(N - 1), 0)
mat_L = np.array([np.insert(upperDiag, 0, 0), mainDiag, np.insert(lowerDiag, int(N), 0)])

# right matrix
mainDiag = -2. * alpha * np.ones(int(N - 1))
mainDiag += 2 * np.ones(int(N - 1))
mainDiag = np.insert(mainDiag, 0, 1)
mainDiag = np.insert(mainDiag, len(mainDiag), 1)
upperDiag = np.insert(alpha * np.ones(int(N - 1)), 0, 0)
lowerDiag = np.insert(alpha * np.ones(int(N - 1)), int(N - 1), 0)
mat_R = np.diag(mainDiag, 0)
mat_R += np.diag(upperDiag, 1)
mat_R += np.diag(lowerDiag, -1)

# boundary temperature for all times
u[0] = u[-1] = 0.
u[0] = u[-1] = 1.

u = np.array(list(map(lambda x: u, T)))

for k in range(0, int(steps)):
    u[k + 1] = linalg.solve_banded((1, 1), mat_L, mat_R.dot(u[k]))

filtered = int(steps / 100)

# plt.plot(X, u[0])
# plt.plot(X, u[100])
# plt.plot(X, u[200])
# plt.plot(X, u[int(steps)])
# plt.plot(X, np.exp(-D * np.pi ** 2 * t) * np.sin(np.pi * X), color='k')

x, t = np.meshgrid(X, T[0::filtered])
plt.pcolormesh(x, t, u[0::filtered])
plt.colorbar()
plt.tight_layout()
plt.savefig('plots/crank_heat.pdf')
