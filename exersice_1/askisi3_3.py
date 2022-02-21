import numpy as np

def seidel(a, x ,b):
    n = len(a)                   
    for j in range(n):        
        d = b[j]                  
        for i in range(0, n):     
            if(j != i):
                d-=a[j][i] * x[i]
        x[j] = d / a[j][j]
    return x    

def inf_norm(matrix):
    return np.linalg.norm(matrix, np.inf)

print("Gauss-Seidel method\n")
n = 10   # 10x10 matrix                            
A = [[0 for x in range(n)] for y in range(n)]
b = [1 for x in range(n)]
x = [0 for i in range(n)]
for i in range(n):
    A[i][i]=5
for i in range(n-1):
    A[i][i+1]=-2
    A[i+1][i]=-2
b[0]=3
b[n-1]=3

print("Solution for 10x10 (n=10):")
e = [pow(10,-4) for x in range(n)]
i=1
while(inf_norm(e)>0.5*pow(10,-4)):
    x_prev = x.copy()
    x = seidel(A, x, b)
    for j in range (n):
        e[j]=x[j]-x_prev[j]
    i+=1
print(x)  
print("Total iterations for 10x10 (n=10):",i,"\n")

n = 10000   # 10000x10000 matrix                            
A = [[0 for x in range(n)] for y in range(n)]
b = [1 for x in range(n)]
x = [0 for i in range(n)]
for i in range(n):
    A[i][i]=5
for i in range(n-1):
    A[i][i+1]=-2
    A[i+1][i]=-2
b[0]=3
b[n-1]=3

print("Solution for 10000x10000 (n=10000):")
e = [pow(10,-4) for x in range(n)]
i=1
while(inf_norm(e)>0.5*pow(10,-4)):
    x_prev = x.copy()
    x = seidel(A, x, b)
    for j in range (n):
        e[j]=x[j]-x_prev[j]
    i+=1
print(x)  
print("Total iterations for 10000x10000:",i)