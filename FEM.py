import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as la
import scipy.sparse as sp
import concurrent.futures

import time

def laplace_solver(V, T, bndry_V, bndry_vals):
    """
    Solves the Laplace equation with a triangular mesh using the finite element method.

    Parameters:
        V (ndarray): An n x 2 array of node positions.
        T (ndarray): An m x 3 array of triangle vertex indices.
        bndry_V (ndarray): A p x 2 array of boundary node positions.
        bndry_vals (ndarray): A p x 1 array of boundary values.

    Returns:
        An n x 1 array of nodal values.
    """
    # Calculate the area of each triangle
    tri_areas = np.zeros(len(T))
    for i in range(len(T)):
        p1, p2, p3 = V[T[i]]
        a = np.linalg.norm(p2 - p1)
        b = np.linalg.norm(p3 - p2)
        c = np.linalg.norm(p1 - p3)
        s = (a + b + c) / 2
        tri_areas[i] = np.sqrt(s * (s - a) * (s - b) * (s - c))

    # Assemble the stiffness matrix
    n = len(V)
    K = sp.lil_matrix((n, n))
    for i in range(n):
        # Find all triangles that include vertex i
        tris = [j for j in range(len(T)) if i in T[j]]
        for j in tris:
            # Calculate the contribution of triangle j to the local stiffness matrix
            p1, p2, p3 = V[T[j]]
            x1, y1 = p1
            x2, y2 = p2
            x3, y3 = p3
            area = tri_areas[j]
            b = np.array([[y2-y3, y3-y1, y1-y2], [x3-x2, x1-x3, x2-x1]]) / (2 * area)
            local_K = area * np.dot(b.T, b)
            # Add the contribution to the global stiffness matrix
            for k in range(3):
                for l in range(3):
                    K[T[j, k], T[j, l]] += local_K[k, l]


    # Apply boundary conditions to the stiffness matrix and right-hand side
    for i in range(len(bndry_V)):
        K[i, :] = 0
        K[i, i] = 1
    b = np.zeros(n)
    for i in range(len(bndry_V)):
        b[i] = bndry_vals[i]
    # Solve the system of equations
    u = la.spsolve(K.tocsr(), b)
    return u

def laplace_solver_batched(V, T, bndry_V, bndry_vals):
    """
    Solves the Laplace equation with a triangular mesh using the finite element method.

    Parameters:
        V (ndarray): An n x 2 array of node positions.
        T (ndarray): An m x 3 array of triangle vertex indices.
        bndry_V (ndarray): A p x 2 array of boundary node positions.
        bndry_vals (ndarray): A p x 3 array of boundary values for RGB.

    Returns:
        An n x 3 array of nodal values for RGB.
    """
    # Calculate the area of each triangle
    p1, p2, p3 = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    a = np.linalg.norm(p2 - p1, axis=1)
    b = np.linalg.norm(p3 - p2, axis=1)
    c = np.linalg.norm(p1 - p3, axis=1)
    s = (a + b + c) / 2
    tri_areas = np.sqrt(s * (s - a) * (s - b) * (s - c))

    # Assemble the stiffness matrix
    n = len(V)
    I = []
    J = []
    S = []
    for i in range(len(T)):
        # Calculate the contribution of triangle i to the local stiffness matrix
        p1, p2, p3 = V[T[i]]
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        area = tri_areas[i]
        b = np.array([[y2-y3, y3-y1, y1-y2], [x3-x2, x1-x3, x2-x1]]) / (2 * area)
        local_K = area * np.dot(b.T, b)
        # Add the contribution to the global stiffness matrix
        for k in range(3):
            for l in range(3):
                I.append(T[i, k])
                J.append(T[i, l])
                S.append(local_K[k, l])
    
    K = sp.coo_matrix((S, (I, J)), shape=(n, n)).tocsr()

    # Apply boundary conditions to the stiffness matrix and right-hand side
    b = np.zeros((n, bndry_vals.shape[1]))
    for i in range(len(bndry_V)):
        K[i, :] = 0
        K[i, i] = 1
        b[i] = bndry_vals[i]

    # Solve the system of equations
    u = la.spsolve(K, b)
    return u
