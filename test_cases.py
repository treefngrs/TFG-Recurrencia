from src.poly import Poly
from src.exp import Exp
from src.recurrence import Recurrence

print('###########')
nhr4 = Recurrence(1,2,pol=Poly(1,1,0),init=[1,1])
print(nhr4)
print(nhr4.solve())
print(nhr4.generate(10))
print('###########')
r = Recurrence(1,2,init = [1,1], pol=Poly(1,2,1))
print(r)
print(r.solve())
print(r.generate(10))
print('###########')
r = Recurrence(1, init = [1], exp=Exp(2))
print(r)
print(r.solve())
print(r.generate(10))
print('###########')
r = Recurrence(3,-2, init = [1,1], pol=Poly(1,3), exp=Exp(2))
print(r)
print(r.solve())
print(r.generate(10))


