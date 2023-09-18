import sys
import numpy as np
import symnmfC 

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


# Define functions for various goals
def symnmf(k, vectors):
    # Calculate the graph Laplacian W
    n = len(vectors)
    m = len(vectors[0])
    A = sym(vectors,n, m)
    D = ddg(A, n)
    W = norm(A, D, n)
    # Calculate the average of all entries of W
    mean = np.mean(W)

    # Initialize H as described in 1.4.1
    H = np.random.uniform(0, 2 * np.sqrt(mean / k), (len(vectors), k))

    # Call the symnmf() method from the C extension module
    final_H = symnmfC.symnmf(H, W, len(H), k)
    
    # Output the final H matrix
    for row in final_H:
        print(','.join(map(str, row)))
    


def sym(vectors,n, m):
    # Call the sym() method from the C extension module
    similarity_matrix = symnmfC.sym(vectors,n, m)
    print("started eden stuff")
    for row in similarity_matrix:
        formatted_row = [f"{element:.4f}" for element in row]
        print(",".join(formatted_row))




def ddg(vectors, n):
    # Call the ddg() method from the C extension module
    diagonal_degree_matrix = symnmfC.ddg(vectors, n)

    # Output the diagonal degree matrix
    for value in diagonal_degree_matrix:
        print(value)


def norm(vectors):
    # Call the norm() method from the C extension module
    normalized_similarity_matrix = symnmfC.norm(vectors)

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
