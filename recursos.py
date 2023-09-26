import numpy as np
import sympy as sym

x = sym.Symbol('x',real=True)
y = sym.Symbol('y',real=True)

def GetLegendre(n,x,y):
    
    y = (x**2 - 1)**n
    
    poly = sym.diff( y,x,n )/(2**n*np.math.factorial(n))
    
    return poly

def GetLegendreRecursive(n,x):

    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = x
    else:
        poly = ((2*n-1)*x*GetLegendreRecursive(n-1,x)-(n-1)*GetLegendreRecursive(n-2,x))/n
   
    return sym.expand(poly,x)

def GetDLegendre(n,x):
    Pn = GetLegendreRecursive(n,x)
    return sym.diff(Pn,x,1)
def GetNewton(f,df,xn,itmax=10000,precision=1e-14):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
def GetRoots(f,df,x,tolerancia = 10):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
                
    Roots.sort()
    
    return Roots
def GetAllRootsGLeg(n):

    xn = np.linspace(-1,1,100)
    
    Legendre = []
    DLegendre = []
    
    for i in range(n+1):
        Legendre.append(GetLegendreRecursive(i,x))
        DLegendre.append(GetDLegendre(i,x))
    
    poly = sym.lambdify([x],Legendre[n],'numpy')
    Dpoly = sym.lambdify([x],DLegendre[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots

def GetWeightsGLeg(n):

    Roots = GetAllRootsGLeg(n)

    

    DLegendre = []
    
    for i in range(n+1):
        DLegendre.append(GetDLegendre(i,x))
    
    Dpoly = sym.lambdify([x],DLegendre[n],'numpy')
    Weights = 2/((1-Roots**2)*Dpoly(Roots)**2)
    
    return Weights

def integral(f,n):
    raices = GetAllRootsGLeg(n)
    pesos = GetWeightsGLeg(n)

    I = 0
    for i in range(n):
        I += pesos[i]*f(raices[i])
    return I

def integralab(f,n,a,b):
    raices = GetAllRootsGLeg(n)
    pesos = GetWeightsGLeg(n)

    I = 0
    for i in range(n):
        I += pesos[i]*f(raices[i]*(b-a)/2 + (b+a)/2)
    
    I = I*(b-a)/2
    return I

def masa(f,n,a,b,c,d):
    raices = GetAllRootsGLeg(n)
    pesos = GetWeightsGLeg(n)

    I = 0

    for i in range(n):
        for j in range(n):
            I += pesos[i]*pesos[j]*f(raices[i]*(b-a)/2 + (b+a)/2, raices[j]*(d-c)/2 + (c+d)/2)  
    I = I* (b-a)/2 * (d-c)/2
    return I


def cm(f1,n,a,b,c,d):
    raices = GetAllRootsGLeg(n)
    pesos = GetWeightsGLeg(n)
    mas = masa(f1,n,a,b,c,d)
    f = lambda x,y: x*f1(x,y)
    I = 0

    for i in range(n):
        for j in range(n):
            I += pesos[i]*pesos[j]*f(raices[i]*(b-a)/2 + (b+a)/2, raices[j]*(d-c)/2 + (c+d)/2)  
    I = I* (b-a)/2 * (d-c)/2
    return I*1/mas

f = lambda x,y: x+2*y**2 

print(cm(f,5,1,3,1,4))
