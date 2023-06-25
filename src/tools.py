from math import sqrt
import numpy as np

#calcula los factores de un entero n
def factors(n):

	if not n:
		return [0]
	
	if n == 1:
		return [1]

	ret = [1]
	ret2 = [n]
	
	for x in range(2, int(sqrt(abs(n))+1)):
		if not n%x:
			ret.append(x)
			if x**2 != n:
				ret2.append(n//x)
	
	return ret+ret2[::-1]

#calcula el producto de una lista
def prod(l):
	ret = 1
	for x in l:
		ret*=x
	return ret

#calcula el coeficiente binomial n sobre k
def binomial(n, k):
    C = [0]*(k+1)
    C[0] = 1
 
    for i in range(1, n+1):
        j = min(i, k)
        for _ in range(j):
            C[j] += C[j-1]
            j -= 1
 
    return C[k]

#calcula factorial de n
def factorial(n):
	if not n:
		return 1
	return n*factorial(n-1)

#halla las veces que cambia el signo dentro de una lista
def signchanges(xs):
	s = 0
	for i in range(1,len(xs)):
		if xs[i-1]*xs[i]<0:
			s+=1
	return s

#resuelve un sistema lineal con eliminación de Gauss
def gauss(a,b):
	b = b.astype(float)
	n = len(b)
	x = np.zeros(n)
	A = np.hstack((a, np.atleast_2d(b).T))
	
	for i in range(n):
		
		rolls = 0
		while(A[i][i] == 0):
			if rolls>n-i:
				return None 
			A[i:] = np.roll(A[i:], n+1)
			rolls +=1
		
		for j in range(i+1, n):
			r = A[j][i]/A[i][i]
			for k in range(n+1):
				A[j][k] = A[j][k] - r*A[i][k]

	for i in range(n-1,-1,-1):
		x[i] = A[i][n]
	
		for j in range(i+1,n):
			x[i] = x[i] - A[i][j]*x[j]

		x[i] = x[i]/A[i][i]
	
	return x