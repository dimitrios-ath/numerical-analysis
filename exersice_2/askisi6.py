import numpy as np

def f(x):
    return np.sin(x)

def simpson(f, interval, n):
    total = 0
    length = interval[1]-interval[0]
    total += f(x_points[0]) + f(x_points[n])
    sum_1 = 0
    for i in range(1, n // 2):
        sum_1 += f(x_points[2 * i])
    sum_1 *= 2
    total += sum_1
    sum_2 = 0
    for j in range(1, (n // 2) + 1):
        sum_2 += f(x_points[2 * j - 1])
    sum_2 *= 4
    total += sum_2
    return (length / (3 * n)) * total


def simpson_error_bound(interval, n, M):
    return (((interval[1]-interval[0]) ** 5) / (180 * (n ** 4))) * M

def trapezoid(f, interval, N):
    total = 0
    total += f(x_points[0]) + f(x_points[N])
    length = interval[1]-interval[0]
    sum_1 = 0.0
    for i in range(1, N):
        sum_1 += f([x_points[i]])
    sum_1 *= 2
    total += float(sum_1)
    return (length / (2 * N)) * total

def trapezoid_error_bound(interval, N, M):
    return (((interval[1]-interval[0])**3)/(12 * (N ** 2)))*M

actual_val = np.cos(0) - np.cos(np.pi/2)
x_points = np.linspace(0.0, np.pi/2, 11)
interval = [0, np.pi/2]
print("Actual intergration value: " + str(actual_val))
print()
print("Intergration with Simpson method: " + str(simpson(f,interval,10)))
print("Simpson method theoretical error: " + str(simpson_error_bound(interval,10,abs(np.sin(np.pi / 2))))) # max of abs fourth derivative of f=sin(x) in [0,pi/2] is at pi/2
print("Simpson method actual error: " + str(abs(simpson(f,interval,10)-actual_val)))
print()
print("Intergration with Trapezoid method: " + str(trapezoid(f,interval,10)))
print("Trapezoid method theoretical error: " + str(trapezoid_error_bound(interval,10,abs(-np.sin(np.pi / 2))))) # max of abs second derivative of f=sin(x) in [0,pi/2] is at pi/2
print("Trapezoid method actual error: " + str(abs(trapezoid(f,interval,10)-actual_val)))