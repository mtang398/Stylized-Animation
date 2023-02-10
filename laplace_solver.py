import numpy as np

def laplace_solver(nodes, elements, bc, bc_value):
    num_nodes = nodes.shape[0]
    num_elements = elements.shape[0]
    # Create the stiffness matrix
    K = np.zeros((num_nodes, num_nodes))
    # Create the load vector
    f = np.zeros(num_nodes)
    for element in elements:
        x1, y1 = nodes[element[0]]
        x2, y2 = nodes[element[1]]
        x3, y3 = nodes[element[2]]
        
        # Compute the area of the triangle
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        
        # Compute the local stiffness matrix
        k = np.array([[(y2-y3)**2, + (x3 - x2)**2, (y3-y1)(y2-y3)+(x1-x3)(x3-x2), (y2-y3)(y1-y2)+(x3-x2)(x2-x1)],
                      [(y3-y1)(y2-y3)+(x1-x3)(x3-x2),(y3-y1)**2+(x1-x3)**2, (y2-y3)(y1-y2)+(x3-x2)(x2-x1)],
                      [(y2-y3)(y1-y2)+(x3-x2)(x2-x1),(y2-y3)(y1-y2)+(x3-x2)(x2-x1),(y1-y2)**2 + (x2-x1)**2]]) / (2 * area)
        
        # Assemble the global stiffness matrix
        K[element[0], element[0]] += k[0, 0]
        K[element[0], element[1]] += k[0, 1]
        K[element[0], element[2]] += k[0, 2]
        K[element[1], element[0]] += k[1, 0]
        K[element[1], element[1]] += k[1, 1]
        K[element[1], element[2]] += k[1, 2]
        K[element[2], element[0]] += k[2, 0]
        K[element[2], element[1]] += k[2, 1]
        K[element[2], element[2]] += k[2, 2]
        
    # Apply the Dirichlect Boundary Condition
    num_boundary_nodes = bc.shape[0]
    for i in range(num_boundary_nodes):
        index = bc[i]
        value = bc_value[i]