import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

def interpolation(x, y, triangles, z):
    triang = mtri.Triangulation(x, y, triangles)
    # Interpolate to regularly-spaced quad grid.
    xi, yi = np.meshgrid(np.linspace(0, 500, 2000), np.linspace(0, 500, 2000))

    interp_lin = mtri.LinearTriInterpolator(triang, z)
    zi_lin = interp_lin(xi, yi)

    interp_cubic_geom = mtri.CubicTriInterpolator(triang, z, kind='geom')
    zi_cubic_geom = interp_cubic_geom(xi, yi)

    interp_cubic_min_E = mtri.CubicTriInterpolator(triang, z, kind='min_E')
    zi_cubic_min_E = interp_cubic_min_E(xi, yi)

    # Set up the figure
    plt.tricontourf(triang,z)
    plt.triplot(triang, 'ko-')
    plt.colorbar()
    plt.show()