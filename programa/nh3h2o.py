"""Propiedades y equilibrio liquido-vapor de NH3-H2O segun IAPWS G4-01
(Tillner-Roth y Friend, 1998). Motor: paquete 'iapws' (clase H2ONH3) con dos
correcciones aplicadas a la derivada composicional. Validado contra las Tablas
6, 7 y 8 de la guia IAPWS G4-01. Composicion x = fraccion MOLAR de amoniaco.
Presiones en MPa, T en K, h y s en base masica (kJ/kg, kJ/kg-K)."""
import numpy as np
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
    "_phir": [("2*Tc12*1.12455*x**1.12455", "2*Tc12*1.125455*x**1.125455")],
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
    ns = {"exp": exp, "IAPWS95": IAPWS95, "NH3": NH3, "H2ONH3": H2ONH3}
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

def P_of(rho_mol, T, x):
    """Presion [MPa] sin evaluar cp/w (evita NaN dentro de la campana)."""
    M = Mm(x); ph = _o._phir(rho_mol*M, T, x)
    Z = 1 + ph["delta"]*ph["fird"]
    return Z*R_U*T*rho_mol/1000.0

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
            if v > 0: return brentq(g, a, b, xtol=1e-13, rtol=1e-15)
        raise ValueError(f"sin raiz de vapor: T={T:.2f} P={P:.5g} x={x:.4f}")
    b = 70.0
    while not np.isfinite(g(b)) or g(b) < 0:
        b -= 2.0
        if b < 12: raise ValueError(f"sin raiz de liquido: T={T:.2f} P={P:.5g} x={x:.4f}")
    a = b
    for _ in range(300):
        a *= 0.97
        v = g(a)
        if np.isfinite(v) and v < 0: return brentq(g, a, b, xtol=1e-13, rtol=1e-15)
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
