import math
import random
import time

def newton(f,Df,x0,e):
    xn = x0
    n=1
    while(n<=100000):  # maximum of 100000 iterations
        fxn = f(xn)
        Dfxn = Df(xn)
        xn = xn - fxn/Dfxn
        if abs(fxn) < e:
            print("Number of iterations:",n)
            return xn
        n+=1

def secant(f,a,b):
    if f(a)*f(b) < 0:
        a_n = a
        b_n = b
        n=1
        while(n<100000):  # maximum of 100000 iterations
            m_n = a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))
            f_m_n = f(m_n)
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0:  # exact solution found
                return m_n
            if abs(f_m_n)<error:
                print("Number of iterations:",n)
                break
            n+=1
        return a_n - f(a_n)*(b_n - a_n)/(f(b_n) - f(a_n))

def bisection(f,a,b):
    #N = math.ceil((math.log(b-a)-math.log(error))/math.log(2)) # number of iterations needed
    #print("Number of iterations based on form:",N)
    if f(a)*f(b) < 0:
        a_n = a
        b_n = b
        n=1
        while(n<100000):  # maximum of 100000 iterations
            m_n = (a_n + b_n)/2
            f_m_n = f(m_n)
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0: # exact solution found
                print("Number of iterations:",n)
                return m_n
            if abs(f_m_n)<error:
                print("Number of iterations:",n)
                break
            n+=1
        return (a_n + b_n)/2

def newton_modified(f,Df1,Df2,x0,e):
    xn = x0
    n=0
    while(n<=100000):  # maximum of 100000 iterations
        fxn = f(xn)
        if abs(fxn) < e:
            print("Number of iterations:",n)
            return xn
        Df1xn = Df1(xn)
        Df2xn = Df2(xn)
        xn = xn - 1/((Df1xn/fxn)-((1/2)*(Df2xn/Df1xn)))
        n+=1

def secant_modified(f,a,b,xn,xn1,xn2):
    if f(a)*f(b) < 0:
        a_n = a
        b_n = b
        nums = [xn,xn1,xn2]
        n=0
        while(n<100000):  # maximum of 100000 iterations
            x0=nums[0]
            x1=nums[1]
            x2=nums[2]
            q = f(x0)/f(x1)
            r = f(x2)/f(x1)
            s = f(x2)/f(x0)
            x3 = x2 - (r*(r-q)*(x2-x1)+(1-r)*s*(x2-x0))/((q-1)*(r-1)*(s-1))
            f_m_n = f(x3)
            if abs(f_m_n)<error:
                print("Number of iterations:",n)
                break
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = x3
            elif f(b_n)*f_m_n < 0:
                a_n = x3
                b_n = b_n
            elif f_m_n == 0:  # exact solution found
                return x3
            nums.pop(0)
            nums.append(x3)
            n+=1
        return x2 - (r*(r-q)*(x2-x1)+(1-r)*s*(x2-x0))/((q-1)*(r-1)*(s-1))

def bisection_modified_(f,a,b):
    N = math.ceil((math.log(b-a)-math.log(error))/math.log(2)) # number of iterations needed
    random.seed()
    print("Number of iterations:",N)
    rand = random.random()
    if f(a)*f(b) < 0:
        a_n = a
        b_n = b
        for n in range(1,N+1):
            m_n = a_n+((b_n-a_n)*rand)
            f_m_n = f(m_n)
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0: # exact solution found
                return m_n
        return a_n+((b_n-a_n)*rand)

def bisection_modified(f,a,b):
    random.seed()
    if f(a)*f(b) < 0:
        a_n = a
        b_n = b
        n=1
        while(n<100000):  # maximum of 100000 iterations
            rand = random.random()
            m_n = a_n+((b_n-a_n)*rand)
            f_m_n = f(m_n)
            if f(a_n)*f_m_n < 0:
                a_n = a_n
                b_n = m_n
            elif f(b_n)*f_m_n < 0:
                a_n = m_n
                b_n = b_n
            elif f_m_n == 0: # exact solution found
                return m_n
            if abs(f_m_n)<error:
                print("Number of iterations:",n)
                break
            n+=1
        return a_n+((b_n-a_n)*rand)


error = (1/2)*10**(-5)

f = lambda x: 54*x**6+45*x**5-102*x**4-69*x**3+35*x**2+16*x-4
f_derivative1 = lambda x: 324*x**5+225*x**4-408*x**3-207*x**2+70*x+16
f_derivative2 = lambda x: 1620*x**4+900*x**3-1224*x**2-414*x+70

print("\nf(x) = 54x^6+45x^5-102x^4-69x^3+35x^2+16x-4")

print("\nNewton–Raphson modified method:")
print("First root: " + str(round(newton_modified(f,f_derivative1,f_derivative2,-1.5,error),5)))
print("Second root: " + str(round(newton_modified(f,f_derivative1,f_derivative2,-0.5,error),5)))
print("Third root: " + str(round(newton_modified(f,f_derivative1,f_derivative2,0,error),5)))
print("Fourth root: " + str(round(newton_modified(f,f_derivative1,f_derivative2,0.5,error),5)))
print("Fifth root: " + str(round(newton_modified(f,f_derivative1,f_derivative2,1,error),5)))

print("\nBisection modified method:")
print("First root: " + str(round(bisection_modified(f,-2,-1),5)))
print("Second root: " + str(round(bisection_modified(f,0,0.4),5)))
print("Third root: " + str(round(bisection_modified(f,0.4,1),5)))
print("Fourth root: " + str(round(bisection_modified(f,1,2),5)))

print("\nSecant modified method:")
print("First root: " + str(round(secant_modified(f,-2,-1,-1.4,-1.41,-1.5),5)))
print("Second root: " + str(round(secant_modified(f,0,0.4,0,0.2,0.4),5)))
print("Third root: " + str(round(secant_modified(f,0.4,1,0.4,0.45,0.6),5)))
print("Fourth root: " + str(round(secant_modified(f,1,2,1,1.2,2),5)))

print("\nExecuting modified bisection method 10 times for the same root:")
for n in range(0,10):
    bisection_modified(f,-2,-1)

print("\nComparing the time elapsed for a specific root in classic and modified methods:")
start = time.time()
newton(f,f_derivative1,-1.5,error)
end =  time.time()
print("Time elapsed in classic Newton–Raphson method:" , f'{end-start:.20f}')
start = time.time()
newton_modified(f,f_derivative1, f_derivative2, -1.5,error)
end =  time.time()
print("Time elapsed in modified Newton–Raphson method:" , f'{end-start:.20f}',"\n")

start = time.time()
bisection(f,-2,-1)
end =  time.time()
print("Time elapsed in classic bisection method:" , f'{end-start:.20f}')
start = time.time()
bisection_modified(f,-2,-1)
end =  time.time()
print("Time elapsed in modified bisection method:" , f'{end-start:.20f}',"\n")

start = time.time()
secant(f,-2,-1)
end =  time.time()
print("Time elapsed in classic secant method:" , f'{end-start:.20f}')
start = time.time()
secant_modified(f,-2,-1,-1.4,-1.41,-1.5)
end =  time.time()
print("Time elapsed in modified secant method:" , f'{end-start:.20f}')

print()

