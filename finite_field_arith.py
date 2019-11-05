class Field:
    def __init__(self, p, n, f):
        self.size = p**n
        self.irr = f#Polynomial(f)
        self.char = p
        self.dim = n
        self.inverses = inverses(p)
        
    def __repr__(self):
        poly = str(self.irr)
        return "Field of size " + str(self.size) + ", Characteristic " + str(self.char) + ", with irreducible polynomial " + poly

class Polynomial:
    def __init__(self, coefs, char):
        if not any(coefs):
            coefs = [0]
        else:
            for i in range(len(coefs)-1, -1, -1):
                if coefs[i] != 0:
                    coefs = coefs[:i+1]
                    break
        self.coefs = coefs
        if char:
            self.coefs = [coef % char for coef in self.coefs]
        self.dim = len(self.coefs)
        self.char = char
        
    def __repr__(self):
        poly = ""
        for i in range(len(self.coefs) - 1, -1, -1):
            if self.coefs[i] != 0:
                c = ""
                x = ""
                exp = ""
                if self.coefs[i] != 1 or i == 0:
                    c = str(self.coefs[i])
                if i != 0:
                    x = "x"
                if i > 1:
                    exp = str(i)
                    x += "^"
                poly += c+x+exp+"+"
        if poly == "":
            poly = "0"
        else:
            poly = poly[:-1]
        return poly

    def __add__(self, other):
        assert self.char == other.char
        
        length = max(len(self.coefs),len(other.coefs))
        selfcoefs = self.coefs + [0 for i in range(length - len(self.coefs))]
        othercoefs = other.coefs + [0 for i in range(length - len(other.coefs))]
        
        new_coefs = [(selfcoefs[i] + othercoefs[i]) for i in range(length)]

        if self.char:
            new_coefs = [coef % self.char for coef in new_coefs]
            
        return Polynomial(new_coefs, self.char)

    def __neg__(self):
        negs = [-coef for coef in self.coefs]
        if self.char:
            negs = [coef % self.char for coef in negs]
        return Polynomial(negs, self.char)

    def __sub__(self, other):
        return self + (-other)

    def __pow__(self, num):
        assert isinstance(num, int)
        new = Polynomial([1], self.char)
        for i in range(num):
            new = new * self
        return new
    def __truediv__(self, other):
        if isinstance(other, int):
            inv = inverses(self.char)[other]
            return Polynomial([(coef*inv) % self.char for coef in self.coefs],self.char)
        return self.poly_div_mod(other)[0]        
    def __rtruediv__(self, other):
        if isinstance(other, int):
            nom = Polynomial([other], self.char)
            return nom / self
    def __mod__(self, other):
        if isinstance(other, int):
            return Polynomial([coef % other for coef in self.coefs], self.char)
        if other.is_const() and not other.is_zero():
            
            return Polynomial([coef % other.coefs[0] for coef in self.coefs], self.char)
        else:
            return self.poly_div_mod(other)[1]
        
    def poly_div_mod(self, other):
        assert not other.is_zero() and self.char == other.char
        char = self.char
        if self.char:
            invs = inverses(char)
        q = Polynomial([0], char)
        r = Polynomial([coef for coef in self.coefs], char)
        while not r.is_zero() and r.deg() >= other.deg():
            if char:
                t = r.coefs[r.deg()] * invs[other.coefs[other.deg()]]#lead(r)/lead(other)
            else:
                t = r.coefs[r.deg()] / other.coefs[other.deg()]
            
            new_poly_coefs = [0 for i in range(r.deg()+1)]
         #   print(new_poly_coefs, r, other, r.deg())
            new_poly_coefs[r.deg() - other.deg()] = t
            new_poly = Polynomial(new_poly_coefs, char)
            q = q + new_poly
            #print(new_poly.char, other.char, r.char)
            r = r - (new_poly * other)
            if char:
                if q != q % char or r != r % char:
                    print("Check")
                
        return (q,r)
    
    def __eq__(self, other):
        length = max(len(self.coefs),len(other.coefs))
        selfcoefs = self.coefs + [0 for i in range(length - len(self.coefs))]
        othercoefs = other.coefs + [0 for i in range(length - len(other.coefs))]
        
        return selfcoefs == othercoefs

    def is_zero(self):
        return self.coefs == [0 for i in range(len(self.coefs))]
        
    def deg(self):
        if self.is_zero():
            return 0
        
        n = len(self.coefs)
        
        for i in range(n - 1, -1, -1):
            if self.coefs[i] != 0:
                return i
        return 0
    def __rmul__(self, other):
        if isinstance(other, int):
            return self*other
    def is_const(self):
        return len(self.coefs) == 1
    def __mul__(self, other):
        if isinstance(other, int):
            return Polynomial([coef*other for coef in self.coefs], self.char)
        if other.is_const():
            return Polynomial([coef*other.coefs[0] for coef in self.coefs], self.char)
        assert self.char == other.char
        #assert len(self.coefs) == len(other.coefs) and self.field == other.field
        length = max(len(self.coefs),len(other.coefs))
        selfcoefs = self.coefs + [0 for i in range(length - len(self.coefs))]
        othercoefs = other.coefs + [0 for i in range(length - len(other.coefs))]
        
        #n = self.field.dim
        mult_coefs = [0 for i in range(2*length-1)]
        for i in range(length):
            for j in range(length):
                mult_coefs[i+j] += selfcoefs[i]*othercoefs[j]
        if self.char:
            mult_coefs = [coef % self.char for coef in mult_coefs]
            
        return Polynomial(mult_coefs, self.char)

    
    def inv(self, irr):
        assert self.char == irr.char
        char = self.char
        t = Polynomial([0],char)
        newt = Polynomial([1],char)
        r = Polynomial([coef for coef in irr.coefs], char)
        newr = Polynomial([coef for coef in self.coefs], char)
        while not newr.is_zero():
        #    print(r,newr)
            q = r.poly_div_mod(newr)[0]
            r, newr = newr, r - (q*newr)
            t, newt = newt, t - (q*newt)
        if r.deg() > 0:
            print("Either irr is not irreducible or self is multiple of irr")
        
        return t * inverses(char)[r.coefs[0]]
          
    
class Element:
    def __init__(self, poly, field):
        #print(poly.coefs)
        assert len(poly.coefs) <= field.dim
        self.poly = poly#Polynomial([coef % field.char for coef in coefs])
        self.field = field

    def __repr__(self):
        return str(self.poly)
    
    def __add__(self, other):
        assert self.field == other.field
        p = self.field.char
        new_poly = (self.poly + other.poly) % p
        return Element(new_poly, self.field)
    
    def __neg__(self):
        return Element((-self.poly) % self.field.char, self.field)
    
    def __eq__(self, other):
        return self.poly == other.poly and self.field == other.field
    
    def is_zero(self):
        return self.poly.is_zero()

    def deg(self):
        return self.poly.deg()

    def __pow__(self, num):
        return Element((self.poly**num) % self.field.irr, self.field)
    def __mul__(self, other):
        assert self.field == other.field
        return Element((self.poly*other.poly) % self.field.irr, self.field)
    def __truediv__(self, other):
        assert self.field == other.field
        other_inv = other.inv()
        return Element(self.poly * other_inv.poly, self.field)
    def inv(self):
        return Element(self.poly.inv(self.field.irr), self.field)
    #def __div__(self, other):
        
    def mod(self, other):
        #assert not is_zero(other)
        q = 1

def inverses(p):
    dic = {}
    for i in range(p):
        for j in range(p):
            if (i*j) % p == 1:
                dic[i] = j
                break
    return dic
    
f = Polynomial([1, 1, 0, 1], 2)
F = Field(2, 3, f)
a = Element(Polynomial([1, 0, 1], 2), F)
b = Element(Polynomial([1, 0, 0], 2), F)

G = Field(5, 3, f)
c = Element(Polynomial([2, 0, 1], 5), G)

g = Polynomial([1,1,0,1,1], 2)
