
import math
import numpy as np
import matplotlib.pyplot as plt

def apply_to_arrays(func):
    def wrapper(x,*args,**kargs):
        if isinstance(x,np.ndarray) : return np.array([func(i,*args,**kargs) for i in x])
        return func(x,*args,**kargs)
    return wrapper

def natural_spline(points):
    equations = [[1]+[0 for i in range(len(points)-1)],[0 for i in range(len(points)-1)]+[1]]
    solutions = [0,0]
    for i in range(1,len(points)-1):
        h1 = points[i][0] - points[i-1][0]
        h2 = points[i+1][0] - points[i][0]
        f0 = points[i-1][1]
        f1 = points[i][1]
        f2 = points[i+1][1]
        equations.append([h1/6 if M==i-1 else (h1+h2)/3 if M==i else h2/6 if M==i+1 else 0 for M in range(len(points))])
        solutions.append((f2-f1)/h2 - (f1-f0)/h1)
    Mis = np.linalg.solve(np.array(equations),np.array(solutions))
    @apply_to_arrays
    def spline(x):
        nonlocal points,Mis
        for i in range(1,len(points)):
            if x< points[i][0]:
                x0,x1 = points[i-1][0],points[i][0]
                hi = x1-x0
                return Mis[i-1]*(x1-x)**3/(6*hi) + Mis[i]*(x-x0)**3/(6*hi) + (points[i-1][1]-Mis[i-1]*hi**2/6)*(x1-x)/hi + (points[i][1]-Mis[i]*hi**2/6)*(x-x0)/hi
    return spline

def newton_iterpol(points:list((int,int))):
    fxs = []
    xs = []
    coefs = []
    for x,y in points:
        xs.append(x)
        fxs.append(y)
        for i in range(len(fxs)-2,-1,-1):
            fxs[i] = (fxs[i+1]-fxs[i])/(x-xs[i])
        coefs.append(fxs[0])
    fx0 = fxs[0]
    def poli(x):
        nonlocal coefs,xs,fx0
        c = 0
        dx = 1
        for xi,coef in zip(xs,coefs):
            c += dx*coef
            dx *= x-xi
        return c
    return poli





@apply_to_arrays
def F(x):
    return x**2+math.sin(6*x)

points = [(i,F(i)) for i in [-1,-5/7,-3/7,-1/7,1/7,3/7,5/7,1]]
f = natural_spline(points)
f1 = newton_iterpol(points)

x_line = np.linspace(-1,1)
y_vals = f(x_line)
y1_val = f1(x_line)
y2_val = F(x_line)
plt.plot(x_line,y_vals,label="Spline Cúbico")
plt.plot(x_line,y1_val,label="Polinónio interpolador")
plt.plot(x_line,y2_val,label="Função x²+sin(6*x)")
plt.show()

