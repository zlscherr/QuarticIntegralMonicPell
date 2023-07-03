F.<a,b> = QQ[]
F.<a,b> = FractionField(F)
B = a^3*(a-1)*(2*a-1)/(a^2-3*a+1)^2
C = -a*(a-1)*(2*a-1)/(a^2-3*a+1)
E = EllipticCurve([1-C,-B,-B,0,0])
P = E(0,0)
if not (10*P).is_zero():
    raise Exception("bad torsion point")

# make the substitution to change to short Weierstrass form

Ep = E.short_weierstrass_model()
phi = E.isomorphism_to(Ep)
Q = phi(E(0,0))
if not (10*Q).is_zero():
    raise Exception("bad torsion point")

R = Ep.a4()/b^4
S = Ep.a6()/b^6
r = Q[0]/b^2
s = Q[1]/b^3

T.<a,b> = ZZ[]
R = T(R.numerator())/T(R.denominator())
S = T(S.numerator())/T(S.denominator())
l = lcm([i.denominator() for i in r.denominator().coefficients()] + [i.denominator() for i in r.numerator().coefficients()])
r = T(l*r.numerator())/T(l*r.denominator())
l = lcm([i.denominator() for i in s.denominator().coefficients()] + [i.denominator() for i in s.numerator().coefficients()])
s = T(l*s.numerator())/T(l*s.denominator())

r2 = (-6*r)(a,b*6)
r1 = (-8*s)(a,b*6)
r0 = (-4*R - 3*r^2)(a,b*6)

R.<x> = T[]
d = x^4 + r2*x^2 + r1*x + r0

load("../continued_fractions.sage")

# sanity check
L = get_cf(d)
if len(L) != 19:
    raise Exception("bad continued fraction")

lc = get_sol(d,L)[0].coefficients()[-1]

print("The following must be integers")
print(8*r2)
print(8*r1)
print(256*r0)
print(lc)

# The values you get from this are
# 8r2 = (-16*a^6 + 64*a^5 - 32*a^4 - 32*a^3 + 16*a - 4)/(a^4*b^2 - 6*a^3*b^2 + 11*a^2*b^2 - 6*a*b^2 + b^2)
# 8r1 = (64*a^5 - 96*a^4 + 32*a^3)/(a^4*b^3 - 6*a^3*b^3 + 11*a^2*b^3 - 6*a*b^3 + b^3)
# 256r0 = (256*a^12 - 2048*a^11 + 7168*a^10 - 14336*a^9 + 16384*a^8 - 6656*a^7 - 6528*a^6 + 9728*a^5 - 4864*a^4 + 768*a^3 + 256*a^2 - 128*a + 16)/(a^8*b^4 - 12*a^7*b^4 + 58*a^6*b^4 - 144*a^5*b^4 + 195*a^4*b^4 - 144*a^3*b^4 + 58*a^2*b^4 - 12*a*b^4 + b^4)
# lc = (a^30*b^20 - 45*a^29*b^20 + 960*a^28*b^20 - 12915*a^27*b^20 + 122955*a^26*b^20 - 881244*a^25*b^20 + 4939025*a^24*b^20 - 22197825*a^23*b^20 + 81410550*a^22*b^20 - 246685725*a^21*b^20 + 623137515*a^20*b^20 - 1320678450*a^19*b^20 + 2359240975*a^18*b^20 - 3563543025*a^17*b^20 + 4560542550*a^16*b^20 - 4950790605*a^15*b^20 + 4560542550*a^14*b^20 - 3563543025*a^13*b^20 + 2359240975*a^12*b^20 - 1320678450*a^11*b^20 + 623137515*a^10*b^20 - 246685725*a^9*b^20 + 81410550*a^8*b^20 - 22197825*a^7*b^20 + 4939025*a^6*b^20 - 881244*a^5*b^20 + 122955*a^4*b^20 - 12915*a^3*b^20 + 960*a^2*b^20 - 45*a*b^20 + b^20)/(-512*a^38 + 6656*a^37 - 40448*a^36 + 152576*a^35 - 400064*a^34 + 773696*a^33 - 1142624*a^32 + 1316096*a^31 - 1196834*a^30 + 864146*a^29 - 495368*a^28 + 224168*a^27 - 79100*a^26 + 21308*a^25 - 4232*a^24 + 584*a^23 - 50*a^22 + 2*a^21)

# This gives you all possible rational values of a.  The idea is to kill off b
# and then use the lemma from the paper
f = (256*r0)^5*lc
a_set = get_vals(f)

# This gives you all possible rational values of b.  We can consider 8r1 + lc
# clear all denominators and then use the same Lemma as above.
ab_set = set()
for A in a_set:
    for B in get_b_set(8*r1,lc,A):
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
