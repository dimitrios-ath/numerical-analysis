import math

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

error = (1/2)*10**(-5)

f = lambda x: 14*x*math.e**(x-2)-12*math.e**(x-2)-7*x**3+20*x**2-26*x+12
print("\nf(x) = 14xe^(x-2)-12e^(x-2)-7x^3+20x^2-26x+12")
print("\nBisection method:")
print("First root: " + str(round(bisection(f,0,1),5)))
print("Second root: " + str(round(bisection(f,1,3),5)))

f_derivative = lambda x: 14*x*math.e**(x-2)+2*math.e**(x-2)-21*x**2+40*x-26
print("\nNewton–Raphson method:")
print("First root: " + str(round(newton(f,f_derivative,1,error),5)))
print("Second root: " + str(round(newton(f,f_derivative,3,error),5)))

print("\nSecant method:")
print("First root: " + str(round(secant(f,0,1),5)))
print("Second root: " + str(round(secant(f,1,3),5)))

print()