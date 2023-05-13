import numpy as np
import matplotlib.pyplot as plt

N = 200
a = 1
b = 1
nxgrid = 100
nygrid = 100
xmin = -2.2
xmax = 2.2
ymin = -2.2
ymax = 2.2

theta = np.linspace(0, 2 * np.pi, N + 1)
theta = theta[:-1]

x = a * np.cos(theta)
y = b * np.sin(theta)

tx = -a * np.sin(theta)
ty = b * np.cos(theta)

ax = -a * np.cos(theta)
ay = -b * np.sin(theta)
nx = ax / np.sqrt(ax ** 2 + ay ** 2)
ny = ay / np.sqrt(ax ** 2 + ay ** 2)
h = 1 / (N) * 2 * np.pi
w = h * np.sqrt(tx ** 2 + ty ** 2)

f = (x**2 - y**2) / (2*np.pi)

K = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        if i == j:
            kappa = (tx[j] * ax[j] - ty[j] * ay[j]) / (tx[j] ** 2 + ty[j] ** 2) ** 1.5
            kappa = 1
            K[i, j] = -kappa / (4 * np.pi) * w[j]
        else:
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dot = nx[j] * dx + ny[j] * dy
            r = np.sqrt(dx ** 2 + dy ** 2)
            K[i, j] = dot / (2 * np.pi * r ** 2)
            K[i, j] = K[i, j] * w[j]

sigma = np.linalg.solve(-np.eye(N) / 2 + K, f)

xlin = np.linspace(xmin, xmax, nxgrid)
ylin = np.linspace(ymin, ymax, nygrid)
X, Y = np.meshgrid(xlin, ylin)

u = np.zeros((nxgrid, nygrid))
for i in range(nxgrid):
    for j in range(nygrid):
        if (X[i, j] / a) ** 2 + (Y[i, j] / b) ** 2 >= 1:
            u[i, j] = np.nan
        else:
            u[i, j] = 0
        for k in range(N):
            dx = X[i, j] - x[k]
            dy = Y[i, j] - y[k]
            dot = nx[k] * dx + ny[k] * dy
            r = np.sqrt(dx ** 2 + dy ** 2)
            u[i, j] = u[i, j] - dot / (2 * np.pi * r ** 2) * w[k] * sigma[k]
            
            
plt.imshow(u)
plt.colorbar()
plt.show()

U = (X**2 - Y**2) / (2*np.pi)
E = u - U + np.finfo(float).eps
plt.contourf(X, Y, np.log10(np.abs(E)))
plt.colorbar()
plt.show()
