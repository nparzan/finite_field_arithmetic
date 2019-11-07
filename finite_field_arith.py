import random


class Field:
    """Finite Field representation.

    A finite field of prime (or prime power) size.

    Attributes:
        char (int): Characteristic
        dim (int): Dimension
        size (int): Number of elements in field
        irr (Polynomial): Irreducible polynomial (if `dim > 1`)
        inverses (dict): Inverses of constants in the field
        indeterminate (str): Indeterminate symbol for polynomials
    """

    def __init__(self, p, n, f, ind="x"):
        """Constructor:

        Args:
            p (int): characteristic
            n (int): dimension
            f (Polynomial): irreducible polynomial
            ind (str): indeterminate symbol (default: "x")


        Returns:
            A Field:
                * If `n == 1`: The field of `p` elements `Fp`
                * If `n > 1`: The field of polynomials in the ring `Fp[x]`
                mudulu the irreducible polynomial `f`: `Fp[x] % <f>`
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
        """Returns the zero element of the field
        """

        return Polynomial([0], self)

    def one(self):
        """Returns the identity element of the field
        """

        return Polynomial([1], self)

    def __repr__(self):
        """Text representation.

        Args:
            self (Field): field

        Returns:
            A description:
                * If `n == 1`: describe field size
                * If `n > 1`: describe field size, characteristic
                and irreducible polynomial
        """
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

    # TODO: evaluation should be in the field.
    # i.e. take every coef and convert to field element and then multiply
    def __call__(self, val):
        exp = 1
        ret = Element(Polynomial([0], self.field), self.field)
        for coef in self.coefs:
            t = Element(Polynomial([coef * exp], self.field), self.field)
            ret += t
            exp = exp * val
            ret = ret % self.field.char
        return ret

    def inv(self, irr):
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
    def __init__(self, poly, field):
        assert len(poly.coefs) <= field.dim
        self.poly = poly
        self.field = field
        self.coef_field = poly.field

    def __repr__(self):
        return str(self.poly) + " in the " + str(self.field).lower()

    def __add__(self, other):
        if isinstance(other, int):
            pol = Polynomial([other], self.coef_field)
            return self + Element(pol, self.field)
        assert self.field == other.field
        p = self.field.char
        new_poly = (self.poly + other.poly) % p
        return Element(new_poly, self.field)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Element((-self.poly) % self.field.char, self.field)

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
        return Element((self.poly**num) % self.field.irr, self.field)

    def __mul__(self, other):
        if isinstance(other, int):
            pol = Polynomial([other], self.coef_field)
            return self * Element(pol, self.field)
        assert self.field == other.field
        mod = self.field.char
        if self.field.irr:
            mod = self.field.irr
        return Element((self.poly * other.poly) % mod, self.field)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, int):
            pol = Polynomial([other], self.coef_field)
            return self / Element(pol, self.field)
        assert self.field == other.field
        other_inv = other.inv()
        el1 = Element(self.poly, self.field)
        el2 = Element(other_inv.poly, self.field)
        return el1 * el2

    def __rtruediv__(self, other):
        if isinstance(other, int):
            return Element(Polynomial([other], self.field), self.field) / self

    def inv(self):
        return Element(self.poly.inv(self.field.irr), self.field)

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
        cand = Element(Polynomial([1], self.coef_field), self.field) * self
        while True:
            if cand in generated:
                break
            generated.add(cand)
            cand = cand * self
        return generated

    def __mod__(self, other):
#        if isinstance(other, int):
#            other_pol = Polynomial([other], self.coef_field)
#            other_el = Element(other_pol, self.field)
#            return self % other_el

        return Element(self.poly % other, self.field)

    @staticmethod
    def random(field):
        p = field.char
        elems = [i for i in range(p)]
        n = field.dim
        coefs = [random.choice(elems) for i in range(n)]
        return Element(Polynomial(coefs, field), field)

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
