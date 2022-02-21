from turtle import color
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
    y = [0 for i in range(n)] # substitution Ly=b
    for i in range(0,n):
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

def least_squares(x, y, degree):
    a = []
    for row in range(len(x)):
        a.append([])
        for column in range(degree):
            a[row].append(x[row]**column)
    return lu(np.matmul(np.transpose(a), a), np.matmul(np.transpose(a), y))

def least_squares_sin(x, degree):
    return_value = 0.0
    for i in range(degree):
        return_value += least_squares_coef[i] * (x**i)
    return return_value

# BTC-USD price prediction

"""
BTC-USD from finance.yahoo.com
Date	        Open	    High	    Low	        Close*	    Adj Close**	Volume
Feb 13, 2021	47,491.20	48,047.75	46,392.28	47,105.52	47,105.52	70,250,456,155
Feb 12, 2021	47,877.04	48,745.73	46,424.98	47,504.85	47,504.85	76,555,041,196
Feb 11, 2021	44,898.71	48,463.47	44,187.76	47,909.33	47,909.33	81,388,911,810
Feb 10, 2021	46,469.76	47,145.57	43,881.15	44,918.18	44,918.18	87,301,089,896
Feb 09, 2021	46,184.99	48,003.72	45,166.96	46,481.11	46,481.11	91,809,846,886
->  Feb 08, 2021	38,886.83	46,203.93	38,076.32	46,196.46	46,196.46	101,467,222,687
Feb 07, 2021	39,250.19	39,621.84	37,446.15	38,903.44	38,903.44	65,500,641,143
Feb 06, 2021	38,138.39	40,846.55	38,138.39	39,266.01	39,266.01	71,326,033,653
Feb 05, 2021	36,931.55	38,225.91	36,658.76	38,144.31	38,144.31	58,598,066,402
Feb 04, 2021	37,475.11	38,592.18	36,317.50	36,926.07	36,926.07	68,838,074,392
Feb 03, 2021	35,510.82	37,480.19	35,443.98	37,472.09	37,472.09	61,166,818,159
Feb 02, 2021	33,533.20	35,896.88	33,489.22	35,510.29	35,510.29	63,088,585,433
Feb 01, 2021	33,114.58	34,638.21	32,384.23	33,537.18	33,537.18	61,400,400,660
Jan 31, 2021	34,270.88	34,288.33	32,270.18	33,114.36	33,114.36	52,754,542,671
Jan 30, 2021	34,295.93	34,834.71	32,940.19	34,269.52	34,269.52	65,141,828,798
"""

x_actual = [i for i in range(1, 16)]
y_actual = [34269.52, 33114.36, 33537.18, 35510.29, 37472.09, 36926.07, 38144.31, 39266.01, 38903.44, 46196.46, 46481.11, 44918.18, 47909.33, 47504.85, 47105.52]

x_data = x_actual[:10]
y_data = y_actual[:10]

for poly_degree in range (2,5):
    least_squares_coef = least_squares(x_data, y_data, (poly_degree+1))

    y_estimated = np.array([None]*len(x_actual))
    for i in range (0, len(x_actual)):
        y_estimated[i] = least_squares_sin(x_actual[i], (poly_degree+1)) # calculate sine values using Newton method

    plt.figure(figsize=(10, 10))
    plt.title("BTC-USD price prediction (" + str(poly_degree) + " degree polynomial)")
    plt.xticks(np.arange(0,16,1))
    plt.grid(True, which='both')
    plt.xlim([1, 15.3])
    plt.ylim([min(np.minimum(y_actual, y_estimated))-1000, max(np.maximum(y_actual, y_estimated))+1000])
    plt.plot(x_actual, y_estimated, label="Predicted price", color='#1e377f')
    plt.plot(x_actual, y_actual, label="Actual price", color='#ff0000')
    plt.scatter(10, y_actual[9])
    plt.scatter(10, y_estimated[9])
    plt.scatter(11, y_actual[10])
    plt.scatter(11, y_estimated[10])
    plt.scatter(15, y_actual[14])
    plt.scatter(15, y_estimated[14])
    plt.xlabel("Day",fontsize=13)
    plt.ylabel("BTC-USD ($)",fontsize=13)
    plt.axvline(10, linestyle='dashed', color="#646464")
    plt.text(9.7,55000,'08/02',rotation=90)
    plt.axvline(11, linestyle='dashed', color="#646464")
    plt.text(10.7,55000,'09/02',rotation=90)
    plt.axvline(15, linestyle='dashed', color="#646464")
    plt.text(14.7,55000,'13/02',rotation=90)                 
    plt.legend()
    plt.show()

# ETH-USD price prediction

"""
ETH-USD from finance.yahoo.com
Date	        Open	    High	    Low	        Close*	    Adj Close**	Volume
Feb 13, 2021	1,843.99	1,871.60	1,770.61	1,814.11	1,814.11	35,359,490,535
Feb 12, 2021	1,783.49	1,861.36	1,744.17	1,843.53	1,843.53	37,905,036,865
Feb 11, 2021	1,743.71	1,806.54	1,708.68	1,783.80	1,783.80	36,021,495,262
Feb 10, 2021	1,768.04	1,826.70	1,686.54	1,744.24	1,744.24	41,916,084,617
Feb 09, 2021	1,746.93	1,815.96	1,711.62	1,768.04	1,768.04	44,180,727,529
->  Feb 08, 2021	1,613.64	1,770.59	1,571.58	1,746.62	1,746.62	48,012,285,956
Feb 07, 2021	1,677.61	1,690.04	1,501.75	1,614.23	1,614.23	39,889,440,151
Feb 06, 2021	1,717.80	1,738.31	1,649.07	1,677.85	1,677.85	39,873,420,648
Feb 05, 2021	1,594.79	1,756.51	1,594.79	1,718.65	1,718.65	40,108,628,454
Feb 04, 2021	1,661.17	1,689.19	1,561.85	1,594.76	1,594.76	44,396,871,836
Feb 03, 2021	1,514.77	1,660.91	1,510.01	1,660.91	1,660.91	41,874,566,399
Feb 02, 2021	1,369.51	1,542.99	1,362.77	1,515.19	1,515.19	45,437,142,801
Feb 01, 2021	1,314.86	1,373.85	1,274.36	1,369.04	1,369.04	29,210,670,920
Jan 31, 2021	1,376.82	1,378.92	1,288.50	1,314.99	1,314.99	25,198,853,581
Jan 30, 2021	1,382.23	1,402.40	1,328.53	1,376.12	1,376.12	30,616,574,234
"""

x_actual = [i for i in range(1, 16)]
y_actual = [1376.12, 1314.99, 1369.04, 1515.19, 1660.91, 1594.76, 1718.65, 1677.85, 1614.23, 1746.62, 1768.04, 1744.24, 1783.80, 1843.53, 1814.11]

x_data = x_actual[:10]
y_data = y_actual[:10]

for poly_degree in range (2,5):
    least_squares_coef = least_squares(x_data, y_data, (poly_degree+1))

    y_estimated = np.array([None]*len(x_actual))
    for i in range (0, len(x_actual)):
        y_estimated[i] = least_squares_sin(x_actual[i], (poly_degree+1)) # calculate sine values using Newton method
    
    plt.figure(figsize=(10, 10))
    plt.title("ETH-USD price prediction (" + str(poly_degree) + " degree polynomial)")
    plt.xticks(np.arange(0,16,1))
    plt.grid(True, which='both')
    plt.xlim([1, 15.3])
    plt.ylim([min(np.minimum(y_actual, y_estimated))-1000, max(np.maximum(y_actual, y_estimated))+1000])
    plt.plot(x_actual, y_estimated, label="Predicted price", color='#1e377f')
    plt.plot(x_actual, y_actual, label="Actual price", color='#ff0000')
    plt.scatter(10, y_actual[9])
    plt.scatter(10, y_estimated[9])
    plt.scatter(11, y_actual[10])
    plt.scatter(11, y_estimated[10])
    plt.scatter(15, y_actual[14])
    plt.scatter(15, y_estimated[14])
    plt.xlabel("Day",fontsize=13)
    plt.ylabel("ETH-USD ($)",fontsize=13)
    plt.axvline(10, linestyle='dashed', color="#646464")
    plt.text(9.7,2500,'08/02',rotation=90)
    plt.axvline(11, linestyle='dashed', color="#646464")
    plt.text(10.7,2500,'09/02',rotation=90)
    plt.axvline(15, linestyle='dashed', color="#646464")
    plt.text(14.7,2500,'13/02',rotation=90)                 
    plt.legend()
    plt.show()