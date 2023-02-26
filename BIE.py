import numpy as np
import matplotlib.pyplot as plt
import time

def solve_laplace_equation(boundary_nodes, boundary_values):
    t1 = time.time()
    n = len(boundary_nodes)

    # Compute the distances between all pairs of boundary nodes
    r = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            r[i, j] = np.sqrt((boundary_nodes[i][0] - boundary_nodes[j][0])**2 + (boundary_nodes[i][1] - boundary_nodes[j][1])**2)

    # Compute the matrix A
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i, j] = np.pi
            else:
                A[i, j] = np.log(r[i, j])

    # Compute the right-hand side vector b
    b = boundary_values - np.mean(boundary_values)

    # Solve the linear system of equations Ax = b
    x = np.linalg.solve(A, b)
    t2 = time.time()
    print('time for solving the linear system of equation')
    print(t2 - t1)
    # Define the solution function
    def solution(x_coord, y_coord):
        u = np.zeros_like(x_coord).astype('float64')
        for i in range(n):
            r = np.sqrt((x_coord - boundary_nodes[i][0])**2 + (y_coord - boundary_nodes[i][1])**2)
            u += x[i] * np.log(r)
        return u + np.mean(boundary_values)

    return solution