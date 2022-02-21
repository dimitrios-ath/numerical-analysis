import math
from pprint import pprint

def cholesky(Α):
    n = len(A)
    L = [[0 for x in range(n)] for y in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum = 0;
            if (j == i):    # summation for diagonals
                for k in range(j):
                    sum += pow(L[j][k], 2);
                L[j][j] = math.sqrt(Α[j][j] - sum);
            else:
                for k in range(j):
                    sum += (L[i][k]*L[j][k]);
                if(L[j][j] > 0):
                    L[i][j] = (Α[i][j] - sum) / L[j][j];
    return L

A = [[4, 12, -16],[12, 37, -43],[-16, -43, 98]];
L = cholesky(A)
print ("Cholesky Decomposition")
print ("A:")
pprint(A)
print ("L:")
pprint(L)
print()