import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def laplace_BIE(node_pos, boundary_condition, eps=1e-8, max_iter=1000):
    """
    Solve Laplace equation using boundary integral equation.

    Parameters:
    -----------
    node_pos : numpy.ndarray of shape (n, 2)
        Array of node positions (x, y) in the boundary.

    boundary_condition : callable
        Function that returns the boundary condition at a given position (x, y).

    eps : float, optional
        Tolerance for stopping criterion. The algorithm stops when the maximum
        difference between two consecutive solutions is less than `eps`.
        
    max_iter : int, optional
        Maximum number of iterations.

    Returns:
    --------
    u : numpy.ndarray of shape (n,)
        Array of the solution at the boundary nodes.
    """
    n = len(node_pos)
    u = np.zeros(n)
    u_old = np.copy(u)

    for k in range(max_iter):
        # Compute the boundary integral
        integral = np.sum([(u[j] - u[i]) / np.linalg.norm(node_pos[i] - node_pos[j]) * boundary_condition(i) 
                           for i in range(n) for j in range(n) if i != j])

        # Compute the new solution
        u = 1 / n * np.sum([boundary_condition(node_pos[i]) for i in range(n)]) + 1 / (2 * np.pi) * integral

        # Check the stopping criterion
        if np.max(np.abs(u - u_old)) < eps:
            break
        else:
            u_old = np.copy(u)

    return u
