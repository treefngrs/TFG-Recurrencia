from poly import Poly
from exp import Exp
from fractions import Fraction
from recurrence import Recurrence
from tools import factors,binomial,factorial

print()
p1 = Poly(1,0,-9,-4,12)
print(p1)
p2 = Poly(1,-3,-10)

print(p1.findallroots())
print(p2.findallroots())


print()
print()
r1 = Recurrence(1,1, init = [1,1])
r2 = Recurrence(3,10)
r3 = Recurrence(3,10, init = [4,1])
r4 = Recurrence(0,9,4,-12)
#print(r1)
#print(r1.generate(10, [2,3]))
print(r2.generate(10))
print('--------')
print(r3.solve())
print('--------')
print(r4.solve())



print()
print(Fraction(0.5))
#print(factors(18))
print()
print()
p2 = Poly(1,0,-2)
print(p2)
print(p2.findratroot())
p3 = Poly(1,0,-6,0,8)
print(p3)
print(p3.findallroots())

print('------')
print(binomial(10,2))
print('###########')
nhr = Recurrence(2,1,pol=Poly(1,0,0,var='n'))
print(nhr)
print(nhr.solve())
nhr2 = Recurrence(2,1,init=[1,1],pol=Poly(1,0,0,var='n'))
print(nhr2)
print(nhr2.solve())
print('###########')
nhr3 = Recurrence(2,-1,pol=Poly(2,var='n'))
print(nhr3)
print(nhr3.solve())
print('###########')
nhr4 = Recurrence(2,-1,exp=Exp(2),init=[1,1])
print(nhr4)
print(nhr4.solve())
print(nhr4.generate(20))
print('###########')
exp1=Exp(2,5,-3,2)
print(exp1)
print(exp1.eval(5))

