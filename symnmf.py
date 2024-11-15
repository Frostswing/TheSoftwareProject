import sys
import numpy as np
import symnmfC


# Set the random seed for consistent initialization
np.random.seed(0)

# Functions Area of Code
def print_and_exit():
    print("An Error Has Occurred")
    sys.exit(1)

def read_file_to_array(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            vector = [float(num) for num in line.strip().split(',')]
            data.append(vector)
    return data

# Goals Functions Area of Code
def symnmf(k, vectors):
    n = len(vectors)
    m = len(vectors[0])
    A = sym(vectors, n, m)
    D = ddg(A, n)
    W = norm(A, n)
    meanW = np.mean(W)

    # Initialize H as described in 1.4.1
    H_np = np.random.uniform(0, 2 * np.sqrt(meanW / k), (len(vectors), k))
    H = H_np.tolist()

    # Call the symnmf() method from the C extension module
    final_H = symnmfC.symnmf(H, W, len(W), k)

    return final_H

def sym(vectors, n, m):
    # Call the sym() method from the C extension module
    similarity_matrix = symnmfC.sym(vectors, n, m)
    return similarity_matrix

def ddg(A, n):
    # Call the ddg() method from the C extension module
    diagonal_degree_matrix = symnmfC.ddg(A, n)
    return diagonal_degree_matrix

def norm(A, n):
    # Call the norm() method from the C extension module
    normalized_similarity_matrix = symnmfC.norm(A, n)
    return normalized_similarity_matrix

# Main Area of Code
def main():
    # Arguments Area of Code
    if len(sys.argv) != 4:
        print_and_exit()
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]

    vectors = read_file_to_array(file_name)  # array of all our vectors
    n = len(vectors)
    m = len(vectors[0])

    if goal == "symnmf":
        output_matrix = symnmf(k, vectors)
    elif goal == "sym":
        output_matrix = sym(vectors, n, m)
    elif goal == "ddg":
        A = sym(vectors, n, m)
        output_matrix = ddg(A, n)
    elif goal == "norm":
        A = sym(vectors, n, m)
        output_matrix = norm(A, n)
    else:
        print_and_exit()
    
    # Prints the output matrix
    for row in output_matrix:
        formatted_row = [f"{element:.4f}" for element in row]
        print(",".join(formatted_row))

if __name__ == "__main__":
    main()