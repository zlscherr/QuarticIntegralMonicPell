# the continued fraction of sqrt(d(x)), assuming it is known to be periodic.
def get_cf(d):
    F = d.constant_coefficient().parent()
    R.<x> = LaurentSeriesRing(F,70)
    d = R(d(x^-1))
    f = d^(1/2)
    v = f.valuation()
    L = [sum(f[k]*x^k for k in range(v,1))]
    while L[-1] != 2*L[0]:
        f = 1/(f - L[-1])
        L.append(sum(f[k]*x^k for k in range(v,1)))
    for i in range(len(L)):
        L[i] = L[i](x^-1)
    return L

# Returns the minimal degree solution for the polynomial Pell equation.  Here d(x)
# Needs to be Pellian for this to succeed.
def get_sol(d, L=None):
    if L == None:
        L = get_cf(d)
    P = [0,1]
    Q = [1,0]
    for a in L:
        P.append(a*P[-1]+P[-2])
        Q.append(a*Q[-1]+Q[-2])
        if (P[-1]^2 - d*Q[-1]^2 == 1):
            R = d.parent()
            L = P[-1].coefficients()
            f = sum(R(x^i)*L[i] for i in range(len(L)))
            L = Q[-1].coefficients()
            g = sum(R(x^i)*L[i] for i in range(len(L)))
            return (f,g)

# Function to get the possible values of a using the method of Lemma 5.7
# in the paper.  This supposes we have a rational function of the form
# p(x)/q(x) where deg(p) > deg(q) and q(0) = 0.  Then this produces
# all rational a for which p(a)/q(a) is an integer.
def get_vals(f):
    R = f.parent().base_ring()
    R.<xx> = R[]
    g = f.substitute(dict([i,xx] for i in f.parent().gens()))
    p = g.numerator()
    q = g.denominator()
    l = lcm([c.denominator() for c in p.coefficients()] + [c.denominator() for c in q.coefficients()])
    p = p*l
    q = q*l
    if q(0) != 0 or p(0) == 0 or q.degree() >= p.degree():
        raise Exception("Cannot apply lemma from paper")
    value_set = set()
    for r in ZZ(p.constant_coefficient()).divisors():
        for s in ZZ(p.leading_coefficient()).divisors():
            if q(r/s) != 0 and p(r/s)/q(r/s) in ZZ:
                value_set.add(r/s)
            if q(-r/s) != 0 and q(-r/s)/q(-r/s) in ZZ:
                value_set.add(-r/s)
    return value_set

# For a fixed value of a=A coming from the above function, this function attempts to find the corresponding
# values of b which make expressions of the form (r*a^m)*(s*b^n) an integer where r,s are integers and a = A
# is a fixed rational number.
def get_b_set(f,g,A):
    #if len(f.denominator().coefficients()) != 1 or len(g.numerator().coefficients()) != 1:
    #    raise Exception("Bad Input")
    S.<x> = ZZ[]
    f = f(A,x)
    g = g(A,x)
    d = f.denominator().degree()
    e = g.numerator().degree()
    b_set = set({1})
    U = f.factor().unit()
    L = 1/g.factor().unit()
    primes = set([i[0] for i in factor(U)] + [i[0] for i in factor(L)])
    for p in primes:
        new = set()
        for exp in range(ceil(1/e*valuation(L,p)),floor(1/d*valuation(U,p))+1):
            for b in b_set:
                new.add(b*p^exp)
        b_set = b_set.union(new)
    return b_set.union({-b for b in b_set})



