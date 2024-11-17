import math
def QR_decomposition(matrix):
    n = len(matrix)
    Q = [[0.0] * n for _ in range(n)]
    R = [[0.0] * n for _ in range(n)]
    for j in range(n):
        q = [matrix[i][j] for i in range(n)]
        for k in range(j):
            R[k][j] = sum(Q[i][k] * matrix[i][j] for i in range(n))
            for i in range(n):
                q[i] -= R[k][j] * Q[i][k]
        norm = math.sqrt(sum(q[i] ** 2 for i in range(n)))
        R[j][j] = norm
        for i in range(n):
            Q[i][j] = q[i] / norm 
    return Q, R

def Hessenberg_reduction(matrix):
    n = len(matrix)
    H = [row[:] for row in matrix] 
    for k in range(n - 2):
        x = [H[i][k] for i in range(k + 1, n)]
        norm_x = math.sqrt(sum(val ** 2 for val in x))
        if norm_x == 0:
            continue
        alpha = -norm_x if x[0] > 0 else norm_x
        v = [x[0] + alpha] + x[1:]
        norm_v = math.sqrt(sum(val ** 2 for val in v))
        v = [val / norm_v for val in v]
        for i in range(k + 1, n):
            for j in range(k, n):
                H[i][j] -= 2 * v[i - (k + 1)] * sum(v[m - (k + 1)] * H[m][j] for m in range(k + 1, n))
        for i in range(n):
            for j in range(k + 1, n):
                H[i][j] -= 2 * v[j - (k + 1)] * sum(H[i][m] * v[m - (k + 1)] for m in range(k + 1, n))
    return H

def QR_algorithm(matrix, max_iterations=100, tolerance=1e-10):
    n = len(matrix) 
    A = [row[:] for row in matrix]  
    for _ in range(max_iterations):
        Q, R = QR_decomposition(A)
        A = [[sum(R[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        off_diagonal = sum(A[i][j] ** 2 for i in range(n) for j in range(i))
        if off_diagonal < tolerance:
            break
    return [A[i][i] for i in range(n)]

def main():
    print("Enter the size of the matrix (n x n):")
    n = int(input())
    if n < 1:
        print("Matrix size must be >=1.")
        return
    print(f"Enter the n x n} matrix row by row, space separated:")
    matrix = []
    for _ in range(n):
        row = list(map(float, input().split()))
        if len(row) != n:
            print("Invalid input. Each row must have exactly n elements.")
            return
        matrix.append(row)

    if n < 50:
        eigenvalues = QR_algorithm(matrix)
    else:
        hessenberg_matrix = essenberg_reduction(matrix)
        eigenvalues = QR_algorithm(hessenberg_matrix)

    print("Eigenvalues are:")
    for eig in eigenvalues:
        print(eig)

main()
