import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect

def lu(A,B): # from first project
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
    y = [0 for i in range(n)] 
    for i in range(0,n): # substitution Ly=b
        y[i] = B[i]/float(L[i][i])
        for k in range(0,i,1):
            y[i] -= y[k]*L[i][k]
    x = [0 for i in range(n)]
    for i in range(n-1,-1,-1):
        x[i] = y[i]
        for j in range(i+1,n):
           x[i] = x[i] - U[i][j]*x[j]
        x[i] = x[i]/U[i][i]
    return x

def newton_coef(x, y):
    n = len(y)
    c = np.zeros([n, n])
    c[:,0] = y
    for j in range(1,n): # divided differences table
        for i in range(n-j):
            c[i][j] = (c[i+1][j-1] - c[i][j-1]) / (x[i+j]-x[i])
    return c

def newton_sin(x):
    n = len(x_data) - 1  # polynomial degree 
    p = newton_coef[n]
    for k in range(1, n+1):
        p = newton_coef[n-k] + (x - x_data[n-k]) * p
    return p

def spline_find_differences(x): # find differences between points of x
    array = [None]*(len(x)-1)
    for i in range (0, len(x) - 1):
        array[i] = x[i+1]-x[i]
    return array

def spline_create_target(n, h, y):
    return [0] + [6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) / (h[i] + h[i-1]) for i in range(1, n - 1)] + [0]

def spline_create_matrix(n, h):
    A = [h[i] / (h[i] + h[i + 1]) for i in range(n - 2)] + [0]
    B = [2] * n
    C = [0] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(n - 2)]
    return A, B, C

def spline_compute_spline(x, y): # compute spline with a given set of points (x, y)
    n = len(x) # number of points
    h = spline_find_differences(x)
    A, B, C = spline_create_matrix(n, h)
    D = spline_create_target(n, h, y)
    M = spline_solve_system(A, B, C, D)
    coefficients = [[(M[i+1]-M[i])*h[i]*h[i]/6, M[i]*h[i]*h[i]/2, (y[i+1] - y[i] - (M[i+1]+2*M[i])*h[i]*h[i]/6), y[i]] for i in range(n-1)]

    def spline_sin(val):
        i = min(bisect(x, val)-1, n-2)
        z = (val - x[i]) / h[i]
        C = coefficients[i]
        return (((C[0] * z) + C[1]) * z + C[2]) * z + C[3]
    return spline_sin

def spline_solve_system(A, B, C, D):
    c_p = C + [0]
    X = [0] * len(B)
    d_p = [0] * len(B)
    d_p[0] = D[0] / B[0]
    c_p[0] = C[0] / B[0]
    for i in range(1, len(B)):
        c_p[i] = c_p[i] / (B[i] - c_p[i - 1] * A[i - 1])
        d_p[i] = (D[i] - d_p[i - 1] * A[i - 1]) / (B[i] - c_p[i - 1] * A[i - 1])
    X[-1] = d_p[-1]
    for i in range(len(B) - 2, -1, -1):
        X[i] = d_p[i] - c_p[i] * X[i + 1]
    return X

def least_squares(x, y):
    degree = len(x)
    a = []
    for row in range(len(x)):
        a.append([])
        for column in range(degree):
            a[row].append(x[row]**column)
    return lu(np.matmul(np.transpose(a), a), np.matmul(np.transpose(a), y))

def least_squares_sin(x):
    return_value = 0.0
    for i in range(len(least_squares_coef)):
        return_value += least_squares_coef[i] * (x ** i)
    return return_value


x_data = np.array(np.linspace(-np.pi, np.pi, 10))
y_data = np.array([None]*len(x_data))
for i in range (0, len(x_data)):
    y_data[i] = np.sin(x_data[i]) # calculate actual sine values for 10 different points in [-pi,pi]

newton_coef = newton_coef(x_data, y_data)[0, :] # calculate coefficients of poly
least_squares_coef = least_squares(x_data, y_data)
spline_sin = spline_compute_spline(x_data, y_data)

x = np.linspace(-np.pi, np.pi, 200)
y_actual = np.array([None]*len(x))
for i in range (0, len(x)):
    y_actual[i] = np.sin(x[i]) # calculate actual sine values for 200 different points in range [-pi, pi]

y_newton = newton_sin(x) # calculate sine values using Newton method

y_splines = np.array([None]*len(x))
for i in range (0, len(x)):
    y_splines[i] = spline_sin(x[i]) # calculate sine values using Spline method

y_least_squares = np.array([None]*len(x))
for i in range (0, len(x)):
    y_least_squares[i] = least_squares_sin(x[i]) # calculate sine values using Newton method

errors_newton = [None]*200
for i in range (0,200):
    errors_newton[i] = abs(y_actual[i]-y_newton[i]) # calculate the interpolation error for newton method

errors_splines = [None]*200
for i in range (0,200):
    errors_splines[i] = abs(y_actual[i]-y_splines[i]) # calculate the interpolation error for splines method

errors_least_squares = [None]*200
for i in range (0,200):
    errors_least_squares[i] = abs(y_actual[i]-y_least_squares[i]) # calculate the interpolation error for least squares method

print("Newton method errors:")
print("Min error:" + str(min(errors_newton)))
print("Max error:" + str(max(errors_newton)))
print("Newton RMSE: " + str(np.sqrt(sum((y_actual - y_newton) ** 2) / len(x))))
print()
print("Spline method errors:")
print("Min error:" + str(min(errors_splines)))
print("Max error:" + str(max(errors_splines)))
print("Splines RMSE: " + str(np.sqrt(sum((y_actual - y_splines) ** 2) / len(x))))
print()
print("Least squares method errors:")
print("Min error:" + str(min(errors_least_squares)))
print("Max error:" + str(max(errors_least_squares)))
print("Least squares RMSE: " + str(np.sqrt(sum((y_actual - y_least_squares) ** 2) / len(x))))


# newton method
plt.figure(figsize=(10, 10))
plt.grid(True, which='both')
plt.axvline(x=0, color='k')
plt.axhline(y=0, color='k')
plt.xlim(-3.5, 3.5)
plt.ylim([-1.1, 1.1])
plt.xlabel("x")
plt.ylabel("sin(x)")
plt.title("Newton polynomial interpolation")
plt.plot(x, y_actual, label='Sine')
plt.plot(x, y_newton, label='Newton interpolation')
plt.legend(loc=2)
plt.show()

# splines method
plt.figure(figsize=(10, 10))
plt.grid(True, which='both')
plt.axvline(x=0, color='k')
plt.axhline(y=0, color='k')
plt.xlim(-3.5, 3.5)
plt.ylim([-1.1, 1.1])
plt.xlabel("x")
plt.ylabel("sin(x)")
plt.title("Spline polynomial interpolation")
plt.plot(x, y_actual, label='Sine')
plt.plot(x, y_splines, label='Splines interpolation')
plt.legend(loc=2)
plt.show()

# least squares method
plt.figure(figsize=(10, 10))
plt.grid(True, which='both')
plt.axvline(x=0, color='k')
plt.axhline(y=0, color='k')
plt.xlim(-3.5, 3.5)
plt.ylim([-1.1, 1.1])
plt.xlabel("x")
plt.ylabel("sin(x)")
plt.title("Least squares polynomial interpolation")
plt.plot(x, y_actual, label='Sine')
plt.plot(x, y_least_squares, label='Least Squares interpolation')
plt.legend(loc=2)
plt.show()

# error plot plot
plt.figure(figsize=(10, 10))
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.xlabel("x")
plt.ylabel("error")
plt.xlim([-3.5, 3.5])
plt.ylim([-0.00001, 0.0007])
plt.title("Interpolation error plot")
plt.plot(x, errors_newton, label='Newton error')
plt.plot(x, errors_splines, label='Splines error')
plt.plot(x, errors_least_squares, label='Least Squares error')
plt.legend(loc=2)
plt.show()
