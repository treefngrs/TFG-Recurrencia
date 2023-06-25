from .tools import factors,prod,signchanges
from fractions import Fraction
from collections import defaultdict
from math import sqrt

class Poly:
    
    def __init__(self, *coeffs, var = 'x'):
        self.coeffs = list(coeffs)
        self.degree = len(self.coeffs)-1
        self.var = var

    #método de clase para copiar un polinomio
    @classmethod
    def copy(cls, poly):
        return cls(*poly.coeffs, var=poly.var)

     
    def __repr__(self):
        return "Poly" + str(tuple(self.coeffs))

    #representación del polinomio en LaTeX
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

    #evalua el polinomio para un valor usando el algoritmo de Hörner
    def eval(self, x):
        ret = self.coeffs[0]
        for c in self.coeffs[1:]:
            ret = ret*x+c
        return ret

    #determina si un valor es raiz del polinomio
    def isroot(self, x):
        return self.eval(x)==0
    
    #halla la n-derivada del polinomio
    def diff(self, n=1):
        newCoeffs=[prod(range(self.degree-i-n+1,self.degree-i+1))*c for i,c in enumerate(self.coeffs[:-n])]
        return Poly(*newCoeffs, var=self.var)

    #division euclidea de dos polinomios (devuelve cociente, resto)
    def longdiv(self, pol2):
        coeffs = self.coeffs.copy()
        quo = []
        for _ in range(self.degree-pol2.degree+1):
            if coeffs[0]==0:
                quo.append(0)
                continue
            f = Fraction(coeffs[0], pol2.coeffs[0])
            quo.append(f)
            coeffs = [c-f*pol2.coeffs[i+1] for i,c in enumerate(coeffs[1:1+pol2.degree])]+coeffs[pol2.degree+1:]
        return (Poly(*quo),Poly(*coeffs))

    #devuelve el polinomio multiplicado por un valor constante
    def mult(self, n):
        coeffs = self.coeffs.copy()
        return Poly(*(c*n for c in coeffs))

    #aisla las raices del polinomio en intervalos usando los polinomios de sturm
    def sturm(self):
        pols = [Poly.copy(self),self.diff()]
        for i in range(2,self.degree+1):
            pols.append(pols[i-2].longdiv(pols[i-1])[1].mult(-1))
        
        n = self.degree
        a,b,c = self.coeffs[0], self.coeffs[1], self.coeffs[2]

        #cota de laguerre para raices reales
        bounds = sorted((-b/(n*a)-((n-1)/(n*a))*sqrt((b**2)-(2*n/(n-1))*a*c), -b/(n*a)+((n-1)/(n*a))*sqrt((b**2)-(2*n/(n-1))*a*c)))
        bounds = [bounds[0]-1,bounds[1]+1]
        v1 = [x for x in (p.eval(bounds[0]) for p in pols) if x != 0]
        v2 = [x for x in (p.eval(bounds[1]) for p in pols) if x != 0]
        v = signchanges(v1)-signchanges(v2)
        bounds = [(*bounds,v)]

        #aislar raices en intervalos
        while(v>1):
            v=0
            for bs in bounds:
                if bs[2]>1:
                    bounds.remove(bs)
                    mid = (bs[0]+bs[1])/2
                    v1 = list(p.eval(bs[0]) for p in pols)
                    v2 = list(p.eval(mid) for p in pols)
                    v3 = signchanges(v1)-signchanges(v2)
                    v = max(v,v3)
                    if v3>0:
                        bounds.append((bs[0],mid,v3))
                    v1 = list(p.eval(bs[1]) for p in pols)
                    v3 = signchanges(v2)-signchanges(v1)
                    v = max(v,v3)
                    if v3>0:
                        bounds.append((mid,bs[1],v3))

        return list((a,b) for a,b,_ in bounds)

    #aproxima por el metodo de newton una raiz a partir de un valor xn
    def newton(self,xn):
        xm = xn - self.eval(xn)/self.diff().eval(xn)
        if abs(xm-xn)<0.0001:
            return xm
        else:
            return self.newton(xm)

    #racionaliza una raiz racional usando Ruffini
    def ruffini(self, r, var2='x'):
        ret = [self.coeffs[0]]
        for c in self.coeffs[1:self.degree]:
            ret.append(c+r*ret[-1])
        return Poly(*ret, var=var2)

    #encuentra una raiz racional usando el teorema de la raiz racional
    def findratroot(self):
        for a in factors(self.coeffs[-1]):
            for b in factors(self.coeffs[0]):
                if self.isroot(Fraction(a,b)):
                    return Fraction(a,b)
                if self.isroot(-1*Fraction(a,b)):
                    return -1*Fraction(a,b)
        return None

    #halla todas las raices del polinomio
    #halla primero las racionales y sus multiplicidades, luego el resto
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
                bounds = pol.sturm()
                for a,b in bounds:
                    ret[pol.newton((a+b)/2)]+=1
                return ret
        return ret


