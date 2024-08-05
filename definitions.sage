S.<a,b,c,d,e,t,r,s,w,g> =QQ[]
S=S.fraction_field()

# Models for the family of elliptic curves E_T. We use the following naming convention:
#    E_{C_N}  is EN and E30 corresponds to T=C_{3}^0


def E3(a,b):
    return EllipticCurve([a,0,a^2*b,0,0])
def E30(a): # This is E_{C_3^0}
    return EllipticCurve([0,0,a,0,0])
def E5(a,b):
    return EllipticCurve([a-b,-a*b,-a^2*b,0,0])
def E7(a,b):
    return EllipticCurve([a^2 + a*b - b^2, a^2*b^2 - a*b^3, a^4*b^2 - a^3*b^3, 0, 0])

# Models for the family of elliptic curves \widetilde{E}_T. We use the following naming convention:
#    \widetilde{E}_{C_N}  is EtN and Et30 corresponds to \widetilde{E}_{C_{3}^0}

def E3t(a,b):
    return EllipticCurve([a,0,a^2*b,-5*a^3*b,-a^4*b*(a+7*b)])
def E30t(a): # This is E_{C_3^0}
    return EllipticCurve([0,0,a,0,-7*a^2])
def E5t(a,b):
    return EllipticCurve([a-b,-a*b,-a^2*b,5*a*b*(a^2-2*a*b-b^2),a*b*(a^4-15*a^3*b+5*a^2*b^2-10*a*b^3-b^4)])
def E7t(a,b):
    return EllipticCurve([a^2+a*b-b^2,a*b^2*(a-b),a^3*b^2*(a-b),5*b*a^7-35*b^2*a^6+70*b^3*a^5-70*b^4*a^4+35*b^5*a^3-5*b^7*a,b*a^11-19*b^2*a^10+94*b^3*a^9-258*b^4*a^8+393*b^5*a^7-343*b^6*a^6+202*b^7*a^5-107*b^8*a^4+46*b^9*a^3-8*b^10*a^2-b^11*a])

#Models for the family of elliptic curve F_{C_3} and \widetilde{F}_{C_3} in the article. Note that these models are obtained from E_{C_3} and \widetilde{E}_{C_3}, respectively, via the isomorphisms [c^2d,0,0,0].

def F3(c,d,e,b):
    return E3(c^3*d^2*e,b).change_weierstrass_model(c^2*d)
def F3t(c,d,e,b):
    return E3t(c^3*d^2*e,b).change_weierstrass_model(c^2*d)

def Fac(g): ## if g=0, return 0. Otherwise, the output is the factorization of g
    try:
        factor(g)
    except ArithmeticError:
        return 0
    else:
        return factor(g)

def Ais(E): ## returns to the factorizations of the Weierstrass coefficients of E
    print('a1=', Fac(E.a1()))
    print('a2=', Fac(E.a2()))
    print('a3=', Fac(E.a3()))
    print('a4=', Fac(E.a4()))
    print('a6=', Fac(E.a6()))

def sig(E): ## Returns the signature (c4,c6,\Delta) of an elliptic curve E
    return Fac(E.c4()), Fac(E.c6()), Fac(E.discriminant())


