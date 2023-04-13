from tools import factors,prod
from fractions import Fraction
from collections import defaultdict
import numpy as np

class Poly:
    
    def __init__(self, *coeffs, var = 'x'):
        self.coeffs = list(coeffs)
        self.degree = len(self.coeffs)-1
        self.var = var

    @classmethod
    def copy(cls, poly):
        return cls(*poly.coeffs, var=poly.var)

     
    def __repr__(self):
        return "Poly" + str(tuple(self.coeffs))

    
    def __str__(self):

        if not self.coeffs:
            return "0"
        res = ""
        for i,c in enumerate(self.coeffs):

            if c==0:
                continue
            
            if c>0 and i>0:
                res += "+"
            
            res += f"{c}"
            
            if self.degree-i > 0:
                if abs(c)==1:
                    res = res.strip('1')
                res += self.var
            
            if self.degree-i > 1:
                res += f"^{self.degree-i}"

        return res

    #algoritmo de Hörner
    def eval(self, x):
        ret = self.coeffs[0]
        for c in self.coeffs[1:]:
            ret = ret*x+c
        return ret

    def isroot(self, x):
        return self.eval(x)==0

    def ruffini(self, r, var2='x'):
        ret = [self.coeffs[0]]
        for c in self.coeffs[1:self.degree]:
            ret.append(c+r*ret[-1])
        return Poly(*ret, var=var2)

    def findratroot(self):
        for a in factors(self.coeffs[-1]):
            for b in factors(self.coeffs[0]):
                if self.isroot(Fraction(a,b)):
                    return Fraction(a,b)
                if self.isroot(-1*Fraction(a,b)):
                    return -1*Fraction(a,b)
        return None

    def findallroots(self):
        #comentar defaultdict
        ret = defaultdict(int)
        pol = Poly.copy(self)
        while(pol.degree>0):
            r = pol.findratroot()
            if r or r==0:
                ret[r]+=1
                pol = pol.ruffini(r)
            else:
                for x in np.roots(pol.coeffs):
                    ret[x]+=1
                return ret
        return ret

    def diff(self, n=1):
        newCoeffs=[prod(range(self.degree-i-n+1,self.degree-i+1))*c for i,c in enumerate(self.coeffs[:-n])]
        return Poly(*newCoeffs, var=self.var)
