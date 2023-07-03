F.<a,b> = QQ[]
F.<a,b> = FractionField(F)
B = a*(a+1)
C = a
E = EllipticCurve([1-C,-B,-B,0,0])
P = E(0,0)
if not (6*P).is_zero():
    raise Exception("bad torsion point")

# make the substitution to change to short Weierstrass form

Ep = E.short_weierstrass_model()
phi = E.isomorphism_to(Ep)
Q = phi(E(0,0))
if not (6*Q).is_zero():
    raise Exception("bad torsion point")

R = Ep.a4()/b^4
S = Ep.a6()/b^6
r = Q[0]/b^2
s = Q[1]/b^3

T.<a,b> = ZZ[]
R = T(R.numerator())/T(R.denominator())
S = T(S.numerator())/T(S.denominator())
r = T(r.numerator())/T(r.denominator())
s = T(s.numerator())/T(s.denominator())

r2 = (-6*r)(a,b*6)
r1 = (-8*s)(a,b*6)
r0 = (-4*R - 3*r^2)(a,b*6)

R.<x> = T[]
d = x^4 + r2*x^2 + r1*x + r0

load("../continued_fractions.sage")

# sanity check
L = get_cf(d)
if len(L) != 11:
    raise Exception("bad continued fraction")

lc = get_sol(d,L)[0].coefficients()[-1]

print("The following must be integers")
print(8*r2)
print(8*r1)
print(256*r0)
print(lc)

# The values you get from this are
# 8r2 = (12*a^2 + 24*a - 4)/b^2
# 8r1 = (32*a^2 + 32*a)/b^3
# 256r0 = (144*a^4 + 64*a^3 + 480*a^2 + 320*a + 16)/b^4
# lc = b^12/(2*a^9 + 8*a^8 + 12*a^7 + 8*a^6 + 2*a^5)

# This gives you all possible rational values of a.  The idea is to kill off b
# and use the Lemma from the paper.
f = (8*r2)^6*lc
a_set = get_vals(f)

# This gives you all possible rational values of b.  We can consider 8r1 + lc
# clear all denominators and then use the same Lemma as above.
ab_set = set()
for A in a_set:
    for B in get_b_set(8*r2,lc,A):
        if (8*r2)(A,B) in ZZ and (8*r1)(A,B) in ZZ and (256*r0)(A,B) in ZZ and lc(A,B) in ZZ:
            ab_set.add((A,B))

# now we can look at all pairs (a,b) and check with c in [1,2,3,4]
integer_d_set = set()
for (r,s) in ab_set:
    for t in [1,2,3,4]:
        e = d(x+t/4)
        S.<y> = QQ[]
        e = sum(y^i*QQ((e.coefficients()[i])(r,s)) for i in range(5))
        if all(p in ZZ for p in e.coefficients()) and is_squarefree(e):
            integer_d_set.add(e)

if len(integer_d_set) == 0:
    print("No d found in this order")
