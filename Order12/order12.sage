F.<a,b> = QQ[]
F.<a,b> = FractionField(F)
B = a*(2*a-1)*(3*a^2-3*a+1)*(2*a^2-2*a+1)/(a-1)^4
C = -a*(2*a-1)*(3*a^2-3*a+1)/(a-1)^3
E = EllipticCurve([1-C,-B,-B,0,0])
P = E(0,0)
if not (12*P).is_zero():
    raise Exception("bad torsion point")

# make the substitution to change to short Weierstrass form

Ep = E.short_weierstrass_model()
phi = E.isomorphism_to(Ep)
Q = phi(E(0,0))
if not (12*Q).is_zero():
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
if len(L) != 23:
    raise Exception("bad continued fraction")

lc = get_sol(d,L)[0].coefficients()[-1]

print("The following must be integers")
print(8*r2)
print(8*r1)
print(256*r0)
print(lc)

# The values you get from this are
# 8r2 = (48*a^8 - 480*a^7 + 1344*a^6 - 1872*a^5 + 1488*a^4 - 672*a^3 + 144*a^2 - 4)/(a^6*b^2 - 6*a^5*b^2 + 15*a^4*b^2 - 20*a^3*b^2 + 15*a^2*b^2 - 6*a*b^2 + b^2) 
# 8r1 = (384*a^6 - 960*a^5 + 1088*a^4 - 672*a^3 + 224*a^2 - 32*a)/(a^4*b^3 - 4*a^3*b^3 + 6*a^2*b^3 - 4*a*b^3 + b^3)
# 256r0 = (2304*a^16 - 9216*a^15 + 33792*a^14 - 155136*a^13 + 544256*a^12 - 1314816*a^11 + 2270976*a^10 - 2911744*a^9 + 2835840*a^8 - 2120448*a^7 + 1217536*a^6 - 531328*a^5 + 172160*a^4 - 39680*a^3 + 6016*a^2 - 512*a + 16)/(a^12*b^4 - 12*a^11*b^4 + 66*a^10*b^4 - 220*a^9*b^4 + 495*a^8*b^4 - 792*a^7*b^4 + 924*a^6*b^4 - 792*a^5*b^4 + 495*a^4*b^4 - 220*a^3*b^4 + 66*a^2*b^4 - 12*a*b^4 + b^4)
# lc = (a^37*b^24 - 37*a^36*b^24 + 666*a^35*b^24 - 7770*a^34*b^24 + 66045*a^33*b^24 - 435897*a^32*b^24 + 2324784*a^31*b^24 - 10295472*a^30*b^24 + 38608020*a^29*b^24 - 124403620*a^28*b^24 + 348330136*a^27*b^24 - 854992152*a^26*b^24 + 1852482996*a^25*b^24 - 3562467300*a^24*b^24 + 6107086800*a^23*b^24 - 9364199760*a^22*b^24 + 12875774670*a^21*b^24 - 15905368710*a^20*b^24 + 17672631900*a^19*b^24 - 17672631900*a^18*b^24 + 15905368710*a^17*b^24 - 12875774670*a^16*b^24 + 9364199760*a^15*b^24 - 6107086800*a^14*b^24 + 3562467300*a^13*b^24 - 1852482996*a^12*b^24 + 854992152*a^11*b^24 - 348330136*a^10*b^24 + 124403620*a^9*b^24 - 38608020*a^8*b^24 + 10295472*a^7*b^24 - 2324784*a^6*b^24 + 435897*a^5*b^24 - 66045*a^4*b^24 + 7770*a^3*b^24 - 666*a^2*b^24 + 37*a*b^24 - b^24)/(10319560704*a^55 - 227030335488*a^54 + 2468954898432*a^53 - 17679987376128*a^52 + 93700321247232*a^51 - 391668606959616*a^50 + 1343822403649536*a^49 - 3889011121864704*a^48 + 9681747055147008*a^47 - 21042861186723840*a^46 + 40388399086339584*a^45 - 69076518643625472*a^44 + 106040414102708736*a^43 - 146970215814500352*a^42 + 184792555115702784*a^41 - 211610768385701376*a^40 + 221399479949504256*a^39 - 212189301331596288*a^38 + 186669068207408192*a^37 - 150978679529491520*a^36 + 112399412957123904*a^35 - 77083424797278848*a^34 + 48718287391936288*a^33 - 28378546312469344*a^32 + 15231325222807512*a^31 - 7527732744167440*a^30 + 3422478246187238*a^29 - 1429453134505214*a^28 + 547484316928794*a^27 - 191850366097544*a^26 + 61337891796580*a^25 - 17831551686844*a^24 + 4694030957196*a^23 - 1113317628976*a^22 + 236458560002*a^21 - 44637562250*a^20 + 7420292678*a^19 - 1073577760*a^18 + 133165900*a^17 - 13881084*a^16 + 1182916*a^15 - 79160*a^14 + 3902*a^13 - 126*a^12 + 2*a^11)

# This gives you all possible rational values of a.  The idea is to kill off b
# and then use the lemma from the paper
f = (8*r2)^8*lc
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
