from finite_field_arith import *


F = Field(7, 1, None)
f = Polynomial([0, 3, 6, 8, 3],F)
g = Polynomial([3, 3, 2], F)
#f+f

F2 = Field(2,1,None)
irr = Polynomial([1, 1, 0, 1], F2)
G = Field(2, 3, irr)
r = Polynomial([0,1,0],G)

# print(G)
# print(r) 
t = Element(r)
# print(-t)
# print(t+2)
# #t = Element(Polynomial([1], G))
# print(t**2)
# print(r/r)
# print(r.field)
# print(t.poly == r)
# print((r.inv()*r)%r.field.irr)
# print(t.inv())
# print(1/t * t)
# print(Element(Polynomial([1], G)))
field = t.generated_subgroup()
for x in field:
    print(x.is_gen())
field = list(field)
print(field[5].inv()+1)
print(Element(Polynomial([1],field[0].field)))
K = Field(7, 1, None)
R = Field(2,1,None)
a = Element(Polynomial([1],R))
z = Element(Polynomial([0], K))
print("z,a:", z,a)
print("z inv:",z.inv(),z.inv().field.indeterminate)
GT = Field(2,3,None,ind="y")
T = Polynomial([field[5].poly,field[2].poly,field[4].poly],GT)
One = Polynomial.one(GT)
print("One:",One)
print("T:",T)
print("T**2:",T*T)

F31 = Field(31,1,None)
irr31 = Polynomial([1, 12, 7, 1], F31)
G = Field(31, 3, irr31)

x = Element.random(G)
print(x, x.inv(), x*x.inv())
# print(x)
# print(len(x.generated_subgroup()))


#import time
#t0 = time.perf_counter()

#y = Element.draw_generator(G, halt = 100)
#print(y)
#print("Found in:", time.perf_counter() - t0)

'''
f = Polynomial([1, 1, 0, 1], 2)
F = Field(2, 3, f)
a = Element(Polynomial([1, 0, 1], 2), F)
b = Element(Polynomial([1, 0, 0], 2), F)

G = Field(5, 3, f)
c = Element(Polynomial([2, 0, 1], 5), G)

g = Polynomial([1,1,0,1,1], 2)

H = Field(5, 1, Polynomial([1],5))
'''
