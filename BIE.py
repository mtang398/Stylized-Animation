import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad
from scipy.optimize import fsolve
import scipy as sp
import time

def green_func(x, y, eps=1e-8):
    """
    The modified Green's function 
    """
    r = np.linalg.norm(x - y)
    if r < eps:
        return 0
    else:
        return -np.log(r) / (2*np.pi)

import numpy as np

def G_normal(z, y):
    
    dist = np.linalg.norm(z - y)
    
    # Compute the unit outward normal vector at the point y on the boundary
    n = (z - y) / dist
    
    # Compute the normal derivative of the Green's function at point y
    gnormal = -1 / (2 * np.pi * dist) * n[1]
    
    return gnormal

def solve_laplace_single_representation(nodes, f, G):
    t1 = time.time()
    N = nodes.shape[0]
    A = np.zeros((N, N))
    b = f

    # Construct the matrix A and the vector b
    for i in range(N):
        for j in range(N):
            if i == j:
                A[i,j] = 0.005
            else:
                A[i,j] = G(nodes[i], nodes[j])

    # Solve the linear system Ax = b
    phi = np.linalg.solve(A, b)
    print(phi)
            
    def u(points):
        M = len(points)
        u = np.zeros(M)
        for i in range(M):
            u[i] = np.sum(phi * np.array([G(points[i], y) for y in nodes]))
        return u
    t2 = time.time()
    print('time spend in total with single representation BIE: ')
    print(t2 - t1)
    return u

def compute_double_layer(nodes, f, G_normal):
    # Compute the kernel matrix K
    K = np.zeros((len(nodes), len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i == j:
                K[i, j] = 0
            else:
                K[i, j] = -G_normal(nodes[i], nodes[j])
                
    A = 0.5*np.identity(len(nodes)) + K
    phi = np.linalg.solve(A, f)
    
    def U(points):
        M = len(points)
        u = np.zeros(M)
        for i in range(M):
            u[i] = - np.sum(phi * np.array([G_normal(points[i], y) for y in nodes]))
        return u

    return U

def solve_laplace_double_representation(nodes, f, G, Gnormal, points):
    t1 = time.time()
    N = nodes.shape[0]
    M = points.shape[0]
    A = np.zeros((N, N))
    
    # Compute the matrix A using the Green function G and its normal derivative dGdn
    for i in range(N):
        z_i = nodes[i]
        for j in range(N):
            if i == j:
                A[i,j] = 0
            else:
                z_j = nodes[j]
                A[i, j] = 0.5 * G(z_i, z_j) - Gnormal(z_i, z_j) * G(z_j, z_i)

    # Solve the linear system A*phi = f to obtain the values of phi
    phi = np.linalg.solve(A, f)
    
    M = len(points)
    u = np.zeros(M)
    for i in range(M):
        u[i] = -np.sum(phi * np.array([Gnormal(points[i], y) for y in nodes]))
    t2 = time.time()
    print('time spend in total with double representation BIE: ')
    print(t2 - t1)
    return u
      

def compute_boundary_normals(nodes, values):
    # Compute the number of nodes on the boundary
    n_nodes = len(nodes)

    # Allocate space for the boundary normals
    normals = np.zeros((n_nodes, 2))

    # Loop over all nodes on the boundary
    for i in range(n_nodes):
        # Compute the index of the next and previous nodes on the boundary
        i_prev = (i - 1) % n_nodes
        i_next = (i + 1) % n_nodes

        # Compute the tangent vector at this node
        tangent = nodes[i_next] - nodes[i_prev]

        # Compute the length of the tangent vector
        tangent_len = np.linalg.norm(tangent)

        # Compute the unit tangent vector
        tangent_unit = tangent / tangent_len

        # Compute the normal vector
        normal = np.array([-tangent_unit[1], tangent_unit[0]])

        # Compute the value gradient
        grad = (values[i_next] - values[i_prev]) / tangent_len

        # Compute the boundary normal
        normals[i] = normal * grad

    # Return the boundary normals
    return normals



'''

def setup_HUGQ(boundary_nodes, node_values, boundary_normals):
    """
    Sets up the matrices H, U, G, and Q for solving the discrete problem.
    
    Args:
    node_values: List of node values.
    boundary_nodes: List of boundary nodes.
    boundary_normals: List of boundary normals.
    
    Returns:
    Tuple of matrices H, U, G, and Q.
    """
    n = len(node_values)
    H = np.zeros((n, n))
    G = np.zeros((n, n))
    U = np.zeros((n,))
    Q = np.zeros((n,))
    
    for i in range(n):
        # Set up local variables for node i
        ui = node_values[i]
        xi, yi = boundary_nodes[i]
        ni, ti = boundary_normals[i]
        hi = 0.0
        gi = 0.0
        
        # Compute influence coefficients
        for j in range(n):
            # Skip node i
            if i == j:
                continue
                
            # Set up local variables for node j
            uj = node_values[j]
            xj, yj = boundary_nodes[j]
            nj, tj = boundary_normals[j]
            hj = 0.0
            gj = 0.0
            
            # Compute integrals for H and G
            q_star = ti
            w = ni
            ds = np.sqrt((xi - xj)**2 + (yi - yj)**2)
            hj += q_star * ds
            gj += w * ds
            
            # Add contribution to H and G
            H[i][j] = hj
            G[i][j] = gj
            
        # Set diagonal of H based on Eq. 37
        H[i][i] = np.sum(H[i]) + 0.5
        
        # Set up U and Q
        U[i] = ui
        Q[i] = np.sum(G[i])
    
    return H, U, G, Q

def solve_laplace_eqn(boundary_node, boundary_value):

    # Compute the boundary normals
    normals = compute_boundary_normals(boundary_node, boundary_value)

    # Set up the matrices H, U, G, and Q
    H, U, G, Q = setup_HUGQ(boundary_node, boundary_value, normals)

    # Solve the linear system to obtain the coefficients A
    A = spsolve(H, Q - G @ U)

    # Define the solution function
    def solution(x, y):
        n = len(boundary_node)
        u = 0
        for i in range(n):
            xi, yi = boundary_node[i]
            ni, _ = normals[i]
            ds = np.sqrt((xi - x)**2 + (yi - y)**2)
            u += ni * A[i] / (2 * np.pi * ds)
        return u

    return solution

'''


'''   
boundary_node = np.array([[0,0],[0,1],[1,0],[1,1]])
boundary_value = np.array([1,1,1,1])
normal = compute_boundary_normals(boundary_node, boundary_value)
H, U, G, Q = setup_HUGQ(boundary_node, boundary_value, normal)
print(H, U, G, Q)
'''