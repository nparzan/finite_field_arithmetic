from finite_field_arith import *


F = Field(7, 1, None)
f = Polynomial([0, 3, 6, 8, 3],F)
g = Polynomial([3, 3, 2], F)
#f+f

F2 = Field(2,1,None)
F13 = Field(13, 1, None)
Elems13 = [Element(Polynomial([i], F13), F13) for i in range(13)]
E0 = Element(Polynomial([0], F2), F2)
E1 = Element(Polynomial([1], F2), F2)
irr = Polynomial([E1, E1, E0, E1], F2)
F8 = Field(2, 3, irr)
a = Element(Polynomial([E0, E1, E1], F2), F8)
print(a.inv(),a)

f_x = Polynomial([Elems13[4], Elems13[7], Elems13[9]], F13)
#print(Elems13[6].inv())
#print(f_x, f_x(4), Elems13[7] % Elems13[4])
exit()


irr = Polynomial([1, 1, 0, 1], F2)
print("eval", irr(1))
G = Field(2, 3, irr, ind="y")
gener = Element.draw_generator(G)
print("generator:",gener)
f = Polynomial([0,1,0],F2)
f_el = Element(f, G)
print("TEST", f_el + f_el)
field4 = list(f_el.generated_subgroup())
print("el2",field4[2])
print("el3",field4[3])
poly_f = Polynomial([field4[2], field4[3], field4[3]], G)
print(poly_f)
#print(poly_f(field4[2]))
print("f1 int",f(1))
print("f1 elem",f(Element(Polynomial([1], F2), F2)))
exit()
F7 = Field(7, 1, None)
g = Polynomial([2, 3, 3, 1], F7)
k = Element(Polynomial([3], F7), F7)
print("g:",g)

print(g(k))
print(g(3))


r = Polynomial([1],F2)
# print(G)
# print(r) 
t = Element(r,G)
q = Element(Polynomial([1], F2), G)
print("t:",t / t)

print("tinv:", t.inv())

F7 = Field(7, 1, None)
pol3 = Polynomial([3], F7)
el3 = Element(pol3, F7)

print("Final")
print(el3*3)


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

