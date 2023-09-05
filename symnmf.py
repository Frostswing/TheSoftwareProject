# symnmf.py

import sys
import numpy as np
import symnmf  # Import your C extension module

# Set the random seed for consistent initialization
np.random.seed(0)

# Functions Area of Code

# Define a function to read the input file


def read_file_to_array(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            vector = [float(num) for num in line.strip().split(',')]
            data.append(vector)
    return data

# Define a function to calculate the graph Laplacian W


def calculate_graph_laplacian(A):
    # Calculate the degree matrix D
    D = np.diag(np.sum(A, axis=1))

    # Calculate D^(-1/2)
    D_sqrt_inv = np.linalg.inv(np.sqrt(D))

    # Calculate the graph Laplacian W
    W = np.dot(np.dot(D_sqrt_inv, A), D_sqrt_inv)

    return W

# Define functions for various goals


def symnmf(k, vectors):
    # Calculate the graph Laplacian W
    W = calculate_graph_laplacian(vectors)

    # Calculate the average of all entries of W
    m = np.mean(W)

    # Initialize H as described in 1.4.1
    H = np.random.uniform(0, 2 * np.sqrt(m / k), (len(vectors), k))

    # Call the symnmf() method from the C extension module
    final_H = symnmf.symnmf(H, W, vectors)

    # Output the final H matrix
    for row in final_H:
        print(','.join(map(str, row)))


def sym(vectors):
    # Call the sym() method from the C extension module
    similarity_matrix = symnmf.sym(vectors)

    # Output the similarity matrix
    for row in similarity_matrix:
        print(','.join(map(str, row)))


def ddg(vectors):
    # Call the ddg() method from the C extension module
    diagonal_degree_matrix = symnmf.ddg(vectors)

    # Output the diagonal degree matrix
    for value in diagonal_degree_matrix:
        print(value)


def norm(vectors):
    # Call the norm() method from the C extension module
    normalized_similarity_matrix = symnmf.norm(vectors)

    # Output the normalized similarity matrix
    for row in normalized_similarity_matrix:
        print(','.join(map(str, row)))


# Arguments Area of Code
if len(sys.argv) != 4:
    print("An Error Has Occurred")
    sys.exit(1)
k = int(sys.argv[1])
goal = sys.argv[2]
file_name = sys.argv[3]

vectors = read_file_to_array(file_name)  # array of all our vectors

if goal == "symnmf":
    symnmf(k, vectors)
elif goal == "sym":
    sym(vectors)
elif goal == "ddg":
    ddg(vectors)
elif goal == "norm":
    norm(vectors)
else:
    print("An Error Has Occurred")
    sys.exit(1)
