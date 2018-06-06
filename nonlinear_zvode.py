import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

steps = 1e5
T = 1.0
dt = T / steps

L = 1.
N = 100.
dx = L / N

D = .1
A = -1.
B = 1.


def conj(x):
    return np.conjugate(x)


class Dynamics:
    def __init__(self, folder, *args):
        self.scanParam = args
        self.folder = folder
        self.data = []

    def equations(self, t, y, params):
        I = 1.0j
        a, = params

        rslt = [(D / (dx ** 2)) * (y[-1] + y[1] - 2. * y[0])]
        for l in range(1, int(N)):
            x = (D / (dx ** 2)) * (y[l - 1] + y[l + 1] - 2. * y[l])
            # add non-linearity
            x += A * y[l] + B * (y[l] ** 3)
            rslt.append(x)
        rslt.append((D / (dx ** 2)) * (y[int(N - 1)] + y[0] - 2. * y[int(N)]))

        return rslt

    def evolve(self):
        initial = []
        for x in np.linspace(0, L, N + 1):
            # initial.append(np.cos(2 * np.pi * x) ** 2)
            # initial.append(0.)
            initial.append(1. / np.sqrt(B))

        initial = np.array(initial)

        # initial[0] = 1.0
        # initial[-1] = 1.0
        # initial[0] = 0.
        # initial[-1] = 0.
        # initial[10:20] = 1.
        initial[0] = 1. / np.sqrt(B)
        initial[-1] = 1. / np.sqrt(B)

        t0, y0 = 0.0, initial

        sol = ode(self.equations).set_integrator('zvode')  # less accurate: method='bdf'
        sol.set_initial_value(y0, t0).set_f_params([5.])

        data_t = [t0]
        data = [np.array(y0)]

        counter = 0.
        while sol.successful() and sol.t < T - dt:
            data_t.append(sol.t + dt)
            data.append(sol.integrate(sol.t + dt))
            counter += 1
            if counter % 10000 == 0:
                print(str(int((counter / steps) * 100.)) + '%')

        self.data = np.array(data).T
        # with open('data/' + self.folder + 'data.dat', 'w') as f:
        #     pickle.dump(data, f)
        #     f.close()


# Single evolution
dynamics = Dynamics('')
dynamics.evolve()

data_t = np.linspace(0, T, steps + 1)
X = np.linspace(0, L, N + 1)

plt.plot(X, dynamics.data.T[0])
plt.plot(X, dynamics.data.T[100])
plt.plot(X, dynamics.data.T[1000])
plt.plot(X, dynamics.data.T[5000])
plt.plot(X, dynamics.data.T[int(steps)])
# plt.plot(data_t, dynamics.data[32])
plt.ticklabel_format(useOffset=False)
plt.xlim(0, 1)
plt.tight_layout()
plt.savefig('plots/zvode.pdf')
