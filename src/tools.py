from math import sqrt

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

def prod(l):
	ret = 1
	for x in l:
		ret*=x
	return ret

def binomial(n, k):
    C = [0]*(k+1)
    C[0] = 1
 
    for i in range(1, n+1):
        j = min(i, k)
        for _ in range(j):
            C[j] += C[j-1]
            j -= 1
 
    return C[k]

def factorial(n):
	if not n:
		return 1
	return n*factorial(n-1)