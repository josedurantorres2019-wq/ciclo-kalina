"""Propiedades y equilibrio liquido-vapor de NH3-H2O segun IAPWS G4-01
(Tillner-Roth y Friend, 1998). Motor: paquete 'iapws' (clase H2ONH3) con dos
correcciones aplicadas a la derivada composicional. Validado contra las Tablas
6, 7 y 8 de la guia IAPWS G4-01. Composicion x = fraccion MOLAR de amoniaco.
Presiones en MPa, T en K, h y s en base masica (kJ/kg, kJ/kg-K)."""
import numpy as np
from math import exp
from scipy.optimize import brentq, fsolve
from iapws.ammonia import H2ONH3, NH3
from iapws.iapws95 import IAPWS95

M_W, M_A = IAPWS95.M, NH3.M

# ---------------------------------------------------------------------------
# CORRECCIONES AL PAQUETE 'iapws' (v1.5.5), aplicadas en tiempo de importacion
# para que este modulo sea autocontenido y reproducible: NO se modifica la
# instalacion. Ambas afectan solo a la derivada composicional dPhi_r/dx, que
# gobierna las fugacidades y por tanto TODO el equilibrio liquido-vapor. Las
# Tablas 6 (region monofasica) de la guia se reproducen con o sin ellas; las
# Tablas 7 y 8 (burbuja y rocio) solo se reproducen CON ellas.
#
#  P1. Errata numerica en dTn/dx: el exponente alfa de la Ec. (4) de G4-01 vale
#      1.125455, pero el ultimo termino de la derivada estaba escrito con
#      1.12455 (dos veces). Efecto pequeno.
#  P2. Falta la regla del producto en dDeltaPhi_r/dx. La Ec. (8) tiene la forma
#          DeltaPhi_r = f(x) * [ A(tau,d) + x*B(tau,d) + x^2*C(tau,d) ],
#      luego  d/dx = f'(x)*[A + xB + x^2C] + f(x)*[B + 2xC].
#      El paquete calculaba solo el primer sumando. Es el origen de la
#      desviacion de 0.5-4 % en las presiones de burbuja y rocio (y del
#      comentario 'FIXME ... differ by 1%, an error I can't find' del autor).
# ---------------------------------------------------------------------------
import inspect as _inspect, textwrap as _tw

_SUBS = {
    "_phir": [
        ("2*Tc12*1.12455*x**1.12455", "2*Tc12*1.125455*x**1.125455"),
        # Las instancias de los puros son CONSTANTES dentro de _phir (solo se
        # usan para llamar a _phir(tau, delta)); construir IAPWS95() y NH3()
        # en cada llamada es ~159 000 construcciones de MEoS por caso (~7 % del
        # perfil). Se construyen UNA vez en _aplicar_correcciones() y se
        # inyectan en el namespace de la funcion.
        ("    water = IAPWS95()\n", "    water = _WATER\n"),
        ("    ammonia = NH3()\n", "    ammonia = _AMMONIA\n"),
    ],
    "_Dphir": [
        ("    firx = dfx*n*delta**d*tau**t\n",
         "    firx = dfx*n*delta**d*tau**t\n    _B = 0.0\n    _C = 0.0\n"),
        ("        firx += x*dfx*n*delta**d*tau**t*exp(-delta**c)",
         "        firx += x*dfx*n*delta**d*tau**t*exp(-delta**c)\n"
         "        _B += n*delta**d*tau**t*exp(-delta**c)"),
        ("    firx += x**2*dfx*n*delta**d*tau**t*exp(-delta**c)",
         "    firx += x**2*dfx*n*delta**d*tau**t*exp(-delta**c)\n"
         "    _C += n*delta**d*tau**t*exp(-delta**c)"),
        ('    prop["firx"] = firx', '    prop["firx"] = firx + fx*(_B + 2*x*_C)'),
    ],
}

def _aplicar_correcciones():
    from math import exp
    ns = {"exp": exp, "IAPWS95": IAPWS95, "NH3": NH3, "H2ONH3": H2ONH3,
          "_WATER": IAPWS95(), "_AMMONIA": NH3()}
    for meth, subs in _SUBS.items():
        src = _tw.dedent(_inspect.getsource(getattr(H2ONH3, meth)))
        src = src.replace("@staticmethod\n", "")
        for a, b in subs:
            if src.count(a) != 1:
                raise RuntimeError(
                    f"iapws.{meth}: no se pudo aplicar la correccion; el paquete "
                    f"instalado difiere del validado (iapws 1.5.5).")
            src = src.replace(a, b)
        exec(compile(src, f"<parche {meth}>", "exec"), ns)
        setattr(H2ONH3, meth, staticmethod(ns[meth]) if meth == "_Dphir" else ns[meth])

_aplicar_correcciones()

_o = H2ONH3()
R_U = 8.314471

def w2m(w): return (w/M_A)/((w/M_A)+((1-w)/M_W))
def m2w(x): return x*M_A/(x*M_A+(1-x)*M_W)
def Mm(x):  return (1-x)*M_W + x*M_A

def prop(rho_mol, T, x):
    """rho_mol en mol/dm3. Devuelve dict con P[MPa], h,s,u [base masica], phi_i."""
    M = Mm(x)
    r = _o._prop(rho_mol*M, T, x)
    Z = r["P"]*1e3*M/(R_U*T*rho_mol*M)
    r["Z"] = Z
    r["phiA"] = r["fugNH3"]/Z**2      # coef. de fugacidad NH3 (Tabla 4, /Z)
    r["phiW"] = r["fugH2O"]/Z**2      # coef. de fugacidad H2O
    r["rho_mol"] = rho_mol
    return r

# ---------------------------------------------------------------------------
# P_of ESPECIALIZADO
#
# El hot path del motor es P_of (45 396 llamadas/caso): cada una invoca
# H2ONH3._phir COMPLETO, que pide a cada puro MEoS._phir (seis derivadas) y a
# _Dphir (siete derivadas). Para P = Z*R*T*rho solo hace falta
# fird = dPhi_r/ddelta. Estas funciones replican EXACTAMENTE ese unico termino
# (mismos terminos, mismo orden de operaciones, por lo que el valor coincide
# con el de _phir), podando todo lo que P_of no consume: ~90 % menos de
# aritmetica por evaluacion de presion, que es el grueso del tiempo por caso.
# La red golden master no deberia moverse ni en el ruido.
# ---------------------------------------------------------------------------
def _mk_fird_puro(const):
    """Devuelve la funcion fird(tau, delta) de UN componente puro a partir de
    sus constantes (IAPWS-95 o NH3). Equivale a MEoS._phir podado a fird."""
    nr1, d1, t1 = const.get("nr1", []), const.get("d1", []), const.get("t1", [])
    nr2, d2, t2 = const.get("nr2", []), const.get("d2", []), const.get("t2", [])
    g2, c2 = const.get("gamma2", []), const.get("c2", [])
    nr3, d3, t3 = const.get("nr3", []), const.get("d3", []), const.get("t3", [])
    a3 = const.get("alfa3", []); e3 = const.get("epsilon3", [])
    b3 = const.get("beta3", []); g3 = const.get("gamma3", [])
    nr4, a4, b4 = const.get("nr4", []), const.get("a4", []), const.get("b4", [])
    A4, B4, C4, D4 = const.get("A", []), const.get("B", []), const.get("C", []), const.get("D", [])
    bt4 = const.get("beta4", [])

    def fird(tau, delta):
        fird = 0.0
        for n, d, t in zip(nr1, d1, t1):
            fird += n*d*delta**(d-1)*tau**t
        for n, d, g, t, c in zip(nr2, d2, g2, t2, c2):
            fird += n*exp(-g*delta**c)*delta**(d-1)*tau**t*(d-g*c*delta**c)
        for n, d, t, a, e, b, g in zip(nr3, d3, t3, a3, e3, b3, g3):
            p = n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)
            fird += p*(d/delta-2*a*(delta-e))
        for n, a, b, A, B, C, D, bt in zip(nr4, a4, b4, A4, B4, C4, D4, bt4):
            tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
            F = exp(-C*(delta-1)**2-D*(tau-1)**2)
            Fd = -2*C*F*(delta-1)
            Delta = tita**2+B*((delta-1)**2)**a
            if Delta == 0:
                DeltaBd = 0.0
            else:
                Deltad = (delta-1)*(A*tita*2/bt*((delta-1)**2)**(0.5/bt-1)
                                    + 2*B*a*((delta-1)**2)**(a-1))
                DeltaBd = b*Delta**(b-1)*Deltad
            fird += n*(Delta**b*(F+delta*Fd)+DeltaBd*delta*F)
        return fird
    return fird

_WATER_FIRD = _mk_fird_puro(IAPWS95._constants)
_AMMONIA_FIRD = _mk_fird_puro(NH3._constants)

# Coeficientes del termino de mezcla (departure, Eq 8 de G4-01): _Dphir podado
# a fird. Los tres grupos y el termino x^2 final, en el MISMO orden que _Dphir.
_DEP_NR1 = (-1.855822e-2,)
_DEP_D1 = (4,)
_DEP_T1 = (1.5,)
_DEP_NR2 = (5.258010e-2, 3.552874e-10, 5.451379e-6, -5.998546e-13, -3.687808e-6)
_DEP_D2 = (5, 15, 12, 12, 15)
_DEP_T2 = (0.5, 6.5, 1.75, 15, 6)
_DEP_C2 = (1, 1, 1, 1, 2)
_DEP_NR3 = (0.2586192, -1.368072e-8, 1.226146e-2, -7.181443e-2, 9.970849e-2,
            1.0584086e-3, -0.1963687)
_DEP_D3 = (4, 15, 4, 5, 6, 10, 6)
_DEP_T3 = (-1, 4, 3.5, 0, -1, 8, 7.5)
_DEP_C3 = (1, 1, 1, 1, 2, 2, 2)


def _fird_departure(tau, delta, x):
    fx = x*(1-x**0.5248379)
    fird = _DEP_NR1[0]*_DEP_D1[0]*delta**(_DEP_D1[0]-1)*tau**_DEP_T1[0]
    for n, d, t, c in zip(_DEP_NR2, _DEP_D2, _DEP_T2, _DEP_C2):
        fird += n*exp(-delta**c)*delta**(d-1)*tau**t*(d-c*delta**c)
    for n, d, t, c in zip(_DEP_NR3, _DEP_D3, _DEP_T3, _DEP_C3):
        fird += x*n*exp(-delta**c)*delta**(d-1)*tau**t*(d-c*delta**c)
    fird += x**2*(-0.7777897)*exp(-delta**2)*delta**(2-1)*tau**4*(2-2*delta**2)
    return fird*fx


def _delta_fird(rho, T, x):
    """(delta, fird) de la MEZCLA. delta = rho/rhon como en _phir; fird =
    (1-x)*fird_agua + x*fird_amon + fird_mezcla. Mismas expresiones y orden de
    operaciones que H2ONH3._phir para conservar el valor exacto."""
    Tc12 = 0.9648407/2*(IAPWS95.Tc+NH3.Tc)
    Tn = (1-x)**2*IAPWS95.Tc + x**2*NH3.Tc + 2*x*(1-x**1.125455)*Tc12
    b = 0.8978069
    rhoc1m = IAPWS95.rhoc/(IAPWS95.M/1000)
    rhoc2m = NH3.rhoc/(NH3.M/1000)
    M = (1-x)*IAPWS95.M + x*NH3.M
    rhoc12 = 1/(1.2395117/2*(1/rhoc1m + 1/rhoc2m))
    rhonm = 1/((1-x)**2/rhoc1m + x**2/rhoc2m + 2*x*(1-x**b)/rhoc12)
    rhon = rhonm*M/1000
    tau = Tn/T
    delta = rho/rhon
    fird = (1-x)*_WATER_FIRD(tau, delta) + x*_AMMONIA_FIRD(tau, delta) \
           + _fird_departure(tau, delta, x)
    return delta, fird


def P_of(rho_mol, T, x):
    """Presion [MPa] sin evaluar cp/w (evita NaN dentro de la campana).

    Version especializada (Fase 1a): en vez de pedir las 6-7 derivadas del
    Helmholtz residual que calcula iapws _phir, evalua SOLO fird = dPhi_r/ddelta
    -- la unica que necesita Z = 1 + delta*fird -- con el mismo orden de
    operaciones, asi que el valor es identico al del _phir completo."""
    M = Mm(x)
    delta, fird = _delta_fird(rho_mol*M, T, x)
    Z = 1 + delta*fird
    return Z*R_U*T*rho_mol/1000.0

def _newton_salva(g, lo, hi, flo=None, fhi=None):
    """Newton-Raphson salvaguardado sobre la funcion de residuo g, con raiz
    acotada en el bracket abierto [lo, hi] (g en MPa, lo que devuelve P_of-P).

    Derivada en diferencias HACIA ADELANTE con paso relativo por rama (un solo
    g() extra por iteracion). Red de seguridad: cualquier paso invalido (fuera
    del bracket, salto desproporcionado, residuo que no mejora, valor no
    finito) deriva a biseccion por signo, que nunca abandona el bracket;
    agotado el limite combinado de Newton+bisecciones, brentq es el ultimo
    recurso. Devuelve la raiz con precision equivalente a brentq(xtol=1e-13)
    con ~4-6 evaluaciones de g en vez de ~15."""
    if flo is None:
        flo = g(lo)
    if fhi is None:
        fhi = g(hi)
    if abs(flo) <= 1e-10:
        return lo
    if abs(fhi) <= 1e-10:
        return hi
    if not (np.isfinite(flo) and np.isfinite(fhi)) or flo*fhi > 0:
        # Sin cambio de signo (g de un extremo no finito, p.ej. borde de
        # dominio, o raiz degenerada): mismo brentq que el codigo viejo, que
        # ya fallaba o convergia aqui segun el caso.
        return brentq(g, lo, hi, xtol=1e-13, rtol=1e-15)
    if abs(flo) <= abs(fhi):
        x0, f0 = lo, flo
    else:
        x0, f0 = hi, fhi
    for _ in range(50):
        if abs(f0) <= 1e-10:
            return x0
        if not np.isfinite(f0):
            return brentq(g, lo, hi, xtol=1e-13, rtol=1e-15)
        eps = 1e-8*max(abs(x0), 0.1)
        f1 = g(x0+eps)
        if np.isfinite(f1):
            deriv = (f1-f0)/eps
            if deriv != 0.0 and np.isfinite(deriv):
                x1 = x0 - f0/deriv
                if lo < x1 < hi and abs(x1-x0) <= 0.5*(hi-lo):
                    if (abs(x1-x0) <= 1e-12*max(abs(x1), 1.0)
                            and abs(f0) <= 1e-10):
                        return x1
                    g1 = g(x1)
                    if np.isfinite(g1) and abs(g1) < abs(f0):
                        x0, f0 = x1, g1
                        continue
        # Biseccion por signo (red de seguridad).
        mid = 0.5*(lo+hi)
        fm = g(mid)
        if not np.isfinite(fm):
            return brentq(g, lo, hi, xtol=1e-13, rtol=1e-15)
        if abs(fm) <= 1e-10:
            return mid
        if flo*fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if abs(flo) <= abs(fhi):
            x0, f0 = lo, flo
        else:
            x0, f0 = hi, fhi
    return brentq(g, lo, hi, xtol=1e-13, rtol=1e-15)

def rho_TPx(T, P, x, phase):
    """Raiz de densidad molar [mol/dm3] para (T,P,x). phase 'v' o 'l'."""
    g = lambda r: P_of(r, T, x) - P
    if phase == "v":
        a = P*1e3/(R_U*T)*0.98
        b = a
        for _ in range(200):
            b *= 1.06
            v = g(b)
            if not np.isfinite(v): b /= 1.06; break
            if v > 0: return _newton_salva(g, a, b, fhi=v)
        raise ValueError(f"sin raiz de vapor: T={T:.2f} P={P:.5g} x={x:.4f}")
    b = 70.0
    while not np.isfinite(g(b)) or g(b) < 0:
        b -= 2.0
        if b < 12: raise ValueError(f"sin raiz de liquido: T={T:.2f} P={P:.5g} x={x:.4f}")
    a = b
    for _ in range(300):
        a *= 0.97
        v = g(a)
        if np.isfinite(v) and v < 0: return _newton_salva(g, a, b, flo=v)
        if a < 1.0: break
    raise ValueError(f"sin raiz de liquido: T={T:.2f} P={P:.5g} x={x:.4f}")

_psat_cache = {}
def Psat_pure(T):
    """Presiones de saturacion [MPa] de NH3 y H2O puros a T (para la inicializacion
    de Raoult). Fuera del dominio de cada puro se extrapola con Clausius-Clapeyron."""
    k = round(T, 6)
    if k in _psat_cache: return _psat_cache[k]
    pa = NH3(T=min(T, 405.3), x=0.5).P
    if T > 405.3: pa *= np.exp(-2500.*(1/T - 1/405.3))
    pw = IAPWS95(T=max(min(T, 646.9), 273.17), x=0.5).P
    _psat_cache[k] = (pa, pw)
    return pa, pw

def _flashP(T, z, kind, tol=1e-10, nmax=300):
    """Sustitucion sucesiva con inicializacion de Raoult.
    kind='bub': z = x_liquido -> devuelve (P_burbuja, x_vapor, conv)
    kind='dew': z = x_vapor   -> devuelve (P_rocio,  x_liquido, conv)"""
    pa, pw = Psat_pure(T)
    if kind == "bub":
        P = z*pa + (1-z)*pw
        y = z*pa/P
    else:
        P = 1.0/(z/pa + (1-z)/pw)
        y = z*P/pa
    y = float(np.clip(y, 1e-8, 1-1e-8))
    P = float(np.clip(P, 1e-8, 39.0))
    for _ in range(nmax):
        xl, xv = (z, y) if kind == "bub" else (y, z)
        try:
            L = prop(rho_TPx(T, P, xl, "l"), T, xl)
            V = prop(rho_TPx(T, P, xv, "v"), T, xv)
            KA, KW = L["phiA"]/V["phiA"], L["phiW"]/V["phiW"]
        except Exception:
            return P, y, False
        if not (np.isfinite(KA) and np.isfinite(KW) and KA > 0 and KW > 0):
            return P, y, False
        if kind == "bub":
            S = z*KA + (1-z)*KW
            yn = z*KA/S
        else:
            S = z/KA + (1-z)/KW
            yn = (z/KA)/S
        Pn = P*S if kind == "bub" else P/S
        if not (np.isfinite(Pn) and np.isfinite(yn)): return P, y, False
        yn = float(np.clip(yn, 1e-8, 1-1e-8)); Pn = float(np.clip(Pn, 1e-8, 39.0))
        err = abs(Pn-P)/P + abs(yn-y)
        P, y = 0.6*P+0.4*Pn, 0.6*y+0.4*yn
        if err < tol: return P, y, True
    return P, y, False

def _polish(T, z, kind, P, y):
    """Pulido de Newton sobre el sistema exacto (rho_L, rho_V, y). Necesario cerca
    del locus critico, donde la sustitucion sucesiva se estanca."""
    xl, xv = (z, y) if kind == "bub" else (y, z)
    try:
        v0 = [rho_TPx(T, P, xl, "l"), rho_TPx(T, P, xv, "v"), y]
    except Exception:
        return P, y, False
    def F(v):
        rL, rV, yy = v
        yy = min(max(yy, 1e-9), 1-1e-9)
        a, b = (z, yy) if kind == "bub" else (yy, z)
        L = prop(rL, T, a); V = prop(rV, T, b)
        return [(L["P"]-V["P"])*1e4,
                np.log(max(a*L["phiA"]*L["P"], 1e-300))-np.log(max(b*V["phiA"]*V["P"], 1e-300)),
                np.log(max((1-a)*L["phiW"]*L["P"], 1e-300))-np.log(max((1-b)*V["phiW"]*V["P"], 1e-300))]
    sol, _, ier, _ = fsolve(F, v0, full_output=True, xtol=1e-13)
    if ier != 1: return P, y, False
    yy = float(min(max(sol[2], 1e-9), 1-1e-9))
    a = z if kind == "bub" else yy
    return prop(sol[0], T, a)["P"], yy, True

def bubbleP(T, x, **kw):
    P, y, ok = _flashP(T, x, "bub", **kw)
    Pp, yp, okp = _polish(T, x, "bub", P, y)
    return (Pp, yp, True) if okp else (P, y, ok)

def dewP(T, x, **kw):
    P, y, ok = _flashP(T, x, "dew", **kw)
    Pp, yp, okp = _polish(T, x, "dew", P, y)
    return (Pp, yp, True) if okp else (P, y, ok)

def Tsat_pure(P):
    """Temperaturas de saturacion de los componentes puros a P [MPa].
    Cota rigurosa: la mezcla sin azeotropo satisface Tsat_NH3 < T_bub,T_dew < Tsat_H2O."""
    fa = lambda T: NH3(T=T, x=0.5).P - P
    fw = lambda T: IAPWS95(T=T, x=0.5).P - P
    return brentq(fa, 200., 405.3, xtol=1e-9), brentq(fw, 274., 647.0, xtol=1e-9)

def _satT(P, x, kind, lo=None, hi=None):
    Ta, Tw = Tsat_pure(P)
    lo = lo or Ta + 0.05
    hi = hi or Tw - 0.05
    fl = bubbleP if kind == "bub" else dewP
    f = lambda T: (fl(T, x)[0] - P)
    T = brentq(f, lo, hi, xtol=1e-8)
    Pc, y, ok = fl(T, x)
    return T, y, ok

def bubbleT(P, x, **kw): return _satT(P, x, "bub", **kw)
def dewT(P, x, **kw):    return _satT(P, x, "dew", **kw)

def h_sat(T, P, x, kind):
    """Entalpia masica del estado saturado (burbuja o rocio) de composicion x."""
    ph = "l" if kind == "bub" else "v"
    return prop(rho_TPx(T, P, x, ph), T, x)["h"]
