import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import pylab as pl


def plot_mollweid(matrix):
    ra = np.linspace(-math.pi, math.pi, len(matrix))
    dec = np.linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))

    X, Y = np.meshgrid(ra, dec)
    Z = matrix.T
    plt.figure()
    ax = pl.subplot(111, projection='mollweide')
    fig = ax.contourf(X, Y, Z, 100)

    # ax.set_title('---$(24h)$', fontsize=15)  # , fontweight='bold')
    plt.xlabel(r'$\theta$', fontsize=15)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=15)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    pl.show()


def resolufication(M, n):
    init_shape = np.shape(M)
    M2 = np.empty((n * init_shape[0], n * init_shape[1]))
    print(np.shape(M2))
    for i in range(init_shape[0] - 1):
        for j in range(init_shape[1] - 1):
            # for k in range(init_shape[1]-1):
            for k in range(i * n, (i + 1) * n):
                for h in range(j * n, (j + 1) * n):
                    M2[k][h] = M[i][j]
    return M2


N = 40

f = np.random.randint(low=-10, high=10, size=(N, 2 * N))

plot_mollweid(f)

n = 5
d = resolufication(f, n)
print(d)
plot_mollweid(d)
