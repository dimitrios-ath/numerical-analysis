def luDecomposition(A,B):
    print("LU Decomposition method")
    n=len(A)
    L = [[0 for x in range(n)] for y in range(n)]
    U = [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        for k in range(i, n): # U
            sum = 0
            for j in range(i):
                sum += (L[i][j] * U[j][k])
            U[i][k] = A[i][k] - sum

        for k in range(i, n): # L
            if (i == k):
                L[i][i] = 1
            else:
                sum = 0
                for j in range(i):
                    sum += (L[k][j] * U[j][i])
                L[k][i] = (A[k][i] - sum)/U[i][i]
 
    print("\nA Matrix:")
    for i in range(n):
        for j in range(n):
            print(A[i][j], end="\t")
        print("", end="\n")

    print("\nB Vector:")
    for i in range(n):
        print(B[i], end="\t")

    print("\n\nL Triangular:")
    for i in range(n):
        for j in range(n):
            print(L[i][j], end="\t")
        print("", end="\n")
    print("\nU Triangular:")
    for i in range(n):
        # U
        for j in range(n):
            print(U[i][j], end="\t")
        print("")
    # substitutioan Ly=b
    y = [0 for i in range(n)]
    for i in range(0,n):
        y[i] = B[i]/float(L[i][i])
        for k in range(0,i,1):
            y[i] -= y[k]*L[i][k]
    
    print("\ny Vector:")
    for i in range(n):
        print(y[i], end="\t")
    
    x = [0 for i in range(n)]
    for i in range(n-1,-1,-1):
        x[i] = y[i]
        for j in range(i+1,n):
           x[i] = x[i] - U[i][j]*x[j]
        x[i] = x[i]/U[i][i]

    return x
    
A = [[1,4,2], [1,2,3],[2,1,3]]
B = [11,11,13]

x = luDecomposition(A,B)
print("\n\nSolution:")
for i in range(len(x)):
    print("x"+str(i)+" = "+str(x[i]), end="\n")
print()


