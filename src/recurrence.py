from .poly import Poly
from .exp import Exp
from .tools import binomial,factorial
import numpy as np
from fractions import Fraction

class Recurrence:

	def __init__(self, *coeffs, init=False, pol=None, exp=None, var='a'):
		self.coeffs = list(coeffs)
		self.r = len(self.coeffs)
		self.init = init
		self.pol = pol
		self.exp = exp
		self.var = var

	def __repr__(self):
		return f"Recurrence({list(coeffs)}, init={self.init})"

	def __str__(self):
		res = f"{self.var}_n = "
		for i,c in enumerate(self.coeffs):

			if c==0:
				continue

			if c>0 and i>0:
				res += "+"

			res += f"{c}"

			if abs(c)==1:
				res = res.strip('1')
			res += self.var+"_{n-"+str(i+1)+'}'
			
		if self.pol:
			if isinstance(self.pol, Poly):
				if self.pol.coeffs[0]>0:
					res += "+"
				res+=self.pol.__str__()

		if self.exp:
			if isinstance(self.exp, Exp):
				if self.exp.base>0:
					res+="+"
				res+=self.exp.__str__()

		if self.init:
			for i,x in enumerate(self.init):
				res += f'\n{self.var}_{i} = {x}'


		return res

	def generate(self, n, init2=False):
		if not self.init and not init2:
			return "No initial values were given"
		if init2:
			ret = init2
		else:
			ret = self.init

		i = len(ret)
		if i<len(self.coeffs):
			return "Not enough initial values were given"
		
		for j in range(i,n):
			if self.pol and isinstance(self.pol, Poly):
				ret.append(sum(x*ret[j-k-1] for k,x in enumerate(self.coeffs))+self.pol.eval(j))
			if self.exp and isinstance(self.exp, Exp):
				ret.append(sum(x*ret[j-k-1] for k,x in enumerate(self.coeffs))+self.exp.eval(j))
			else:
				ret.append(sum(x*ret[j-k-1] for k,x in enumerate(self.coeffs)))

		return ret


	def nhp_pol(self,e):
		m = self.pol.degree+e
		A=[[0]*(m+1) for _ in range(m+1)]
		for j in range(m+1):
			A[j][j]+=1
			for i in range(1,self.r+1):
				for k in range(j,m+1):
					A[j][k]-=self.coeffs[i-1]*(-1)**(k-j)*binomial(k,j)*i**(k-j)

		c=self.pol.coeffs[::-1]+[0]*e
		print(A)
		print(c)
		for i in range(self.r-e,-1,-1):
			c[i]-=sum(A[i][j]*c[j] for j in range(i+1,self.r+1))
			c[i]=Fraction(c[i],A[i][i+e])
		
		return c[::-1]+[0]*e

	def nhp_exp(self,e):
		exp = self.exp
		chpol = Poly(1, *[-x for x in self.coeffs])
		while(chpol.eval(exp.base)==0):
			chpol = chpol.ruffini(exp.base)
		print("NUEVO METODO",Fraction(exp.a*exp.base**(self.r-e),factorial(e)*chpol.eval(exp.base)))
		if not e:
			return Fraction(exp.a,(exp.base**-exp.s)*(1-sum(Fraction(x,exp.base**(k+1)) for k,x in enumerate(self.coeffs))))
		if e:
			return Fraction(-exp.a,(exp.base**-exp.s)*(sum(Fraction(x*(-k-1)**e,exp.base**(k+1)) for k,x in enumerate(self.coeffs))))


	def solve(self):
		chpol = Poly(1, *[-x for x in self.coeffs])
		roots = chpol.findallroots()
		nhppol = None
		nhpexp = None

		if self.pol and isinstance(self.pol, Poly):
			nhppol = Poly(*self.nhp_pol(roots[1]), var='n')
		if self.exp and isinstance(self.exp, Exp):
			nhpexp = Exp(self.exp.base,a=self.nhp_exp(roots[self.exp.base]),m=roots[self.exp.base])
		
		if not self.init:
			ret = f"{self.var}_n = "
			i=1
			for r in roots:
				for j in range(roots[r]):
					ret+=f"c_{i+j}"
					if j:
						ret+="*n"
						if j>1:
							ret+=f"^{j}"
					if r<0:
						ret+=f"*({r})^n + "
					else:
						ret+=f"*{r}^n + "
				i+=roots[r]
			
			if nhppol:
				ret+= '('+nhppol.__str__()+')'

			if nhpexp:
				ret+= '('+nhpexp.__str__()+')'

			return ret.strip('+ ')

		A = []
		for i in range(len(self.init)):
			sub=[]
			for r in roots:
				for j in range(roots[r]):
					sub.append(float((i**j)*(r**i)))
			A.append(sub)
		
		#sustituir 
		A = np.array(A)
		b = np.array(self.init)
		if nhppol:
			b = np.array([float(x-nhppol.eval(i)) for i,x in enumerate(b)])
		if nhpexp:
			b = np.array([float(x-nhpexp.eval(i)) for i,x in enumerate(b)])
		X = np.linalg.solve(A,b)
		if all(isinstance(r,int) or isinstance(r,Fraction) for r in roots):
			s = [Fraction.from_float(x).limit_denominator(6**9) for x in X]
		else:
			s = X
		ret = f"{self.var}_n = "
		i=1
		for r in roots:
			for j in range(roots[r]):
				ret+=f"{s[i+j-1]}"
				if j:
					ret+="*n"
					if j>1:
						ret+=f"^{j}"
				ret+=f"*{r}^n + "

			i+=roots[r]
		
		if nhppol:
			ret+= '('+nhppol.__str__()+')'
		if nhpexp:
			ret+= '('+nhpexp.__str__()+')'

		return ret.strip('+ ')		




