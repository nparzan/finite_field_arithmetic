import random


class Field:
    def __init__(self, p, n, f, ind="x"):
        """ Construct a field. Input: p, n, f.
            Output: The field of size p**n defined as follows:
                1. If n == 1: The field of size Fp (Zp)
                2. If n > 1: The quotient Fp[x]/<f>, that is,
                   the elements are polynomials of degree < n
                   with coefficients in F where addition is done mod p
                   and multiplication mod f
        """
        self.irr = None
        if n > 1:
            self.irr = f
        self.char = p
        self.dim = n
        self.size = p**n
        self.inverses = inverses(p)
        self.indeterminate = ind

    def zero(self):
        return Polynomial([0], self)

    def one(self):
        return Polynomial([1], self)

    def __repr__(self):
        if self.dim == 1:
            return "Field of size " + str(self.size)
        poly = str(self.irr)
        st = "Field of size " + str(self.size)
        st += ", Characteristic " + str(self.char)
        st += ", with irreducible polynomial " + poly
        return st


class Polynomial:
    def __init__(self, coefs, field):

        if not any(coefs):
            coefs = [0]
        else:
            for i in range(len(coefs) - 1, -1, -1):
                if coefs[i] != 0:
                    coefs = coefs[:i + 1]
                    break
        for coef in coefs:
            pass
            # TODO: assert coef in field

        self.coefs = coefs
        if field.char:
            self.coefs = [coef % field.char for coef in self.coefs]

        self.dim = len(self.coefs) - 1
        self.field = field

    def __repr__(self):
        poly = ""
        for i in range(len(self.coefs) - 1, -1, -1):
            if self.coefs[i] != 0:
                c, x, exp = "", "", ""
                if self.coefs[i] != 1 or i == 0:
                    c = str(self.coefs[i])

                if i != 0:
                    x = self.field.indeterminate
                if (x != "x" and self.coefs[i] != 1 and
                        not isinstance(self.coefs[i], int)):
                    c = "(" + c + ")"
                if i > 1:
                    exp = str(i)
                    x += "^"
                poly += c + x + exp + "+"
        if poly == "":
            poly = "00"
        return poly[:-1]

    def __add__(self, other):
        if isinstance(other, int):
            return self + Polynomial([other], self.field)

        assert self.field.char == other.field.char
        selfcoefs, othercoefs, length = pad_lists(self.coefs, other.coefs)
        new_coefs = [(selfcoefs[i] + othercoefs[i]) for i in range(length)]

        if self.field.char:
            new_coefs = [coef % self.field.char for coef in new_coefs]

        return Polynomial(new_coefs, self.field)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        negs = [-coef for coef in self.coefs]
        if self.field.char:
            negs = [coef % self.field.char for coef in negs]
        return Polynomial(negs, self.field)

    def __sub__(self, other):
        return self + (-other)

    def __pow__(self, num):
        assert isinstance(num, int)
        new = Polynomial([1], self.field)
        for i in range(num):
            new = new * self
        return new

    def __truediv__(self, other):
        if isinstance(other, int):
            inv = inverses(self.field.char)[other]
            return self * Polynomial([inv], self.field)
        return self.poly_div_mod(other)[0]

    def __rtruediv__(self, other):
        if isinstance(other, int):
            nom = Polynomial([other], self.field)
            return nom / self

    def __mod__(self, other):
        if isinstance(other, int):
            return Polynomial([c % other for c in self.coefs], self.field)

        if other.is_const() and not other.is_zero():
            coefs = [coef % other.coefs[0] for coef in self.coefs]
            return Polynomial(coefs, self.field)
        else:
            return self.poly_div_mod(other)[1]

    def poly_div_mod(self, other):
        assert not other.is_zero(), "Divide by zero"
        err = "Can't divide polynomials from fields of different chars."
        assert self.field.char == other.field.char, err

        char = self.field.char
        if char:
            invs = inverses(char)
        q = Polynomial([0], self.field)
        r = Polynomial([coef for coef in self.coefs], self.field)
        while not r.is_zero() and r.deg() >= other.deg():
            if char:
                t = r.coefs[r.deg()] * invs[other.coefs[other.deg()]]
            else:
                t = r.coefs[r.deg()] / other.coefs[other.deg()]

            new_poly_coefs = [0 for i in range(r.deg() + 1)]
            new_poly_coefs[r.deg() - other.deg()] = t
            new_poly = Polynomial(new_poly_coefs, self.field)
            q = q + new_poly
            r = r - (new_poly * other)
            if char:
                if q != q % char or r != r % char:
                    print("Check")

        return (q, r)

    def __eq__(self, other):
        if isinstance(other, int):
            return self.coefs[0] == other and len(self.coefs) == 1

        selfcoefs, othercoefs, _ = pad_lists(self.coefs, other.coefs)

        return selfcoefs == othercoefs and self.field == other.field

    def is_zero(self):
        return self == 0 or self.coefs == [0 for i in range(len(self.coefs))]

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
            return self * other

    def is_const(self):
        return len(self.coefs) == 1

    def __mul__(self, other):
        if isinstance(other, int):
            return self * Polynomial([other], self.field)
        err = "Field characteristics must agree when multiplying polynomials"
        assert self.field.char == other.field.char, err

        selfcoefs, othercoefs, length = pad_lists(self.coefs, other.coefs)

        # Do the actual work
        mult_coefs = [0 for i in range(2 * length - 1)]
        for i in range(length):
            for j in range(length):
                mult_coefs[i + j] += selfcoefs[i] * othercoefs[j]

        # Take mod field characteristic
        if self.field.char:
            mult_coefs = [coef % self.field.char for coef in mult_coefs]

        return Polynomial(mult_coefs, self.field)

    def inv(self):
        irr = self.field.irr
        char = self.field.char
        t = Polynomial([0], self.field)
        newt = Polynomial([1], self.field)
        r = Polynomial([1], self.field)
        if irr:
            r = Polynomial([coef for coef in irr.coefs], self.field)
        newr = Polynomial([coef for coef in self.coefs], self.field)
        while not newr.is_zero():
            q = r.poly_div_mod(newr)[0]
            r, newr = newr, r - (q * newr)
            t, newt = newt, t - (q * newt)
        err = "Either field.irr is not irreducible "
        err += "or polynomial is multiple of field.irr"
        assert r.deg() == 0, err

        return t * inverses(char)[r.coefs[0]]

    @staticmethod
    def one(field):
        return field.one()


class Element():
    def __init__(self, poly):
        assert len(poly.coefs) <= poly.field.dim
        self.poly = poly
        self.field = poly.field

    def __repr__(self):
        return str(self.poly) + " in the " + str(self.field).lower()

    def __add__(self, other):
        if isinstance(other, int):
            return self + Element(Polynomial([other], self.field))
        assert self.field == other.field
        p = self.field.char
        new_poly = (self.poly + other.poly) % p
        return Element(new_poly)

    def __neg__(self):
        return Element((-self.poly) % self.field.char)

    def __sub__(self, other):
        return self + (-other)

    def __eq__(self, other):
        if isinstance(other, int):
            return self.poly == Polynomial([other], self.field)
        return self.poly == other.poly and self.field == other.field

    def is_zero(self):
        return self.poly.is_zero()

    def deg(self):
        return self.poly.deg()

    def __pow__(self, num):
        return Element((self.poly**num) % self.field.irr)

    def __mul__(self, other):
        if isinstance(other, int):
            return self * Element(Polynomial([other], self.field))
        assert self.field == other.field
        return Element((self.poly * other.poly) % self.field.irr)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, int):
            return self / Element(Polynomial([other], self.field))
        assert self.field == other.field
        other_inv = other.inv()
        return Element(self.poly) * Element(other_inv.poly)

    def __rtruediv__(self, other):
        if isinstance(other, int):
            return Element(Polynomial([other], self.field)) / self

    def inv(self):
        return Element(self.poly.inv())

    def __hash__(self):
        return hash((tuple(self.poly.coefs), tuple(self.field.irr.coefs)))

    def is_gen(self, verbose=False):
        generated = self.generated_subgroup()
        if verbose:
            for gen in generated:
                print(gen)
        return len(generated) == self.field.size - 1

    def generated_subgroup(self):
        generated = set()
        cand = Element(Polynomial([1], self.field)) * self
        while True:
            if cand in generated:
                break
            generated.add(cand)
            cand = cand * self
        return generated

    @staticmethod
    def random(field):
        p = field.char
        elems = [i for i in range(p)]
        n = field.dim
        coefs = [random.choice(elems) for i in range(n)]
        return Element(Polynomial(coefs, field))

    @staticmethod
    def draw_generator(field, halt=-1):
        while halt != 0:
            cand = Element.random(field)
            gen = cand.generated_subgroup()
            if len(gen) == field.size - 1:
                return cand
            halt -= 1
        return None


def pad_lists(l1, l2):
    length = max(len(l1), len(l2))
    l1_pad = l1 + [0 for i in range(length - len(l1))]
    l2_pad = l2 + [0 for i in range(length - len(l2))]
    return l1_pad, l2_pad, length


def inverses(p):
    dic = {}
    for i in range(p):
        for j in range(p):
            if (i * j) % p == 1:
                dic[i] = j
                break
    return dic
