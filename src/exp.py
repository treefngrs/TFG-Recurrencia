class Exp:
    
    def __init__(self, base, a=1, s=0, m=0, var = 'n'):
        self.base = base    
        self.a = a
        self.s = s
        self.m = m
        self.var = var

    def __repr__(self):
        return f"Exp({('',a)[bool(a-1)]},{self.base},{('',s)[bool(s)]})" 

    #representacion de la expresión en LaTeX
    def __str__(self):
        ret = ""
        if self.a and self.a!=1:
            ret += f"{self.a}*"
        if self.m:
            ret+=f"{self.var}"
            if self.m>1:
                ret+=f"^{self.m}"
            ret+="*"
        if self.base<0:
            ret+=f"({self.base})^"
        else:
            ret+=f"{self.base}^"
        if self.s:
            if self.s<0:
                ret+=f"{{{self.var}{self.s}}}"
            else:
                ret+=f"{{{self.var}+{self.s}}}"
        else:
            ret+=f"{self.var}"
        return ret

    #evalua la expresion para un valor n
    def eval(self,n):
        return self.a*n**self.m*self.base**(n+self.s)