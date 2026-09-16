"""Nucleo de resolucion de ciclos Kalina KCS-11.

Capa de propiedades sobre nh3h2o.py (IAPWS G4-01) + capa de resolucion.
Composiciones: en la INTERFAZ son MASICAS (w); el motor trabaja en MOLARES (x).
La conversion es explicita en cada frontera; no hay ninguna llamada al motor con
una fraccion masica.

Unidades: T [K], P [MPa], h [kJ/kg], s [kJ/kg-K], caudales normalizados a m_b=1.
"""
import bisect
import time
import numpy as np
from scipy.optimize import brentq, fsolve
import nh3h2o as ng

TOL_T = 1e-6          # tolerancia de inversion en temperatura [K]
DQ = 1e-3             # margen del criterio O5

# ---------------------------------------------------------------------------
# CAPA DE PROPIEDADES
# ---------------------------------------------------------------------------
_puros = {}
_flash = {}           # P -> lista ordenada de (T, xL, xV), usada como arranque tibio


def Tsat_puros(P):
    """(T_sat NH3, T_sat H2O) puros a P. Cotas rigurosas de la campana binaria."""
    k = round(P, 10)
    if k not in _puros:
        _puros[k] = ng.Tsat_pure(P)
    return _puros[k]


def flash_TP(T, P):
    """Equilibrio a (T,P): devuelve (x_L, x_V) MOLARES, o None si a esa (T,P) no
    existe region bifasica para ninguna composicion.

    Resuelve la igualdad de fugacidades de los dos componentes. Es el mismo
    equilibrio que bubbleP/dewP de nh3h2o, pero parametrizado por (T,P), que es
    como lo necesita el ciclo. Arranque tibio con el flash previo mas cercano en
    T a la misma presion: reduce las iteraciones de fsolve de ~9 a ~3.
    """
    Ta, Tw = Tsat_puros(P)
    # Banda de guarda de 0.3 K junto a los puros: alli la ventana bifasica en
    # composicion se cierra (x_L -> 1 o x_V -> 0) y el sistema de fugacidades se
    # vuelve singular. Fuera de la banda el estado es monofasico para cualquier
    # composicion salvo un entorno de medida despreciable, y se resuelve como tal.
    if T <= Ta + 0.3 or T >= Tw - 0.3:
        return None
    kP = round(P, 10)
    tabla = _flash.setdefault(kP, [])
    Ts = [r[0] for r in tabla]
    i = bisect.bisect_left(Ts, T)
    if i < len(Ts) and abs(Ts[i] - T) < 1e-9:
        return tabla[i][1], tabla[i][2]
    pa, pw = ng.Psat_pure(T)
    raoult = (float(np.clip((P - pw) / (pa - pw), 1e-4, 1 - 1e-4)),)
    raoult = (raoult[0], float(np.clip(raoult[0] * pa / P, 1e-4, 1 - 1e-4)))
    cand = [tabla[j] for j in (i - 1, i) if 0 <= j < len(tabla)]
    arranques = []
    if cand:
        _, xL0, xV0 = min(cand, key=lambda r: abs(r[0] - T))
        arranques.append((xL0, xV0))          # arranque tibio
    arranques.append(raoult)                  # respaldo: inicializacion de Raoult

    def F(v):
        a = min(max(v[0], 1e-9), 1 - 1e-9)
        b = min(max(v[1], 1e-9), 1 - 1e-9)
        L = ng.prop(ng.rho_TPx(T, P, a, 'l'), T, a)
        V = ng.prop(ng.rho_TPx(T, P, b, 'v'), T, b)
        return [np.log(a * L['phiA']) - np.log(b * V['phiA']),
                np.log((1 - a) * L['phiW']) - np.log((1 - b) * V['phiW'])]

    sol = None
    for v0 in arranques:
        try:
            s, _, ier, _ = fsolve(F, list(v0), full_output=True, xtol=1e-12)
        except ValueError:
            continue
        if ier == 1 and 0 < s[0] < s[1] < 1:
            sol = s
            break
    if sol is None:
        return None
    xL, xV = float(sol[0]), float(sol[1])
    bisect.insort(tabla, (T, xL, xV))
    return xL, xV


def _mono(T, P, x, fase):
    r = ng.prop(ng.rho_TPx(T, P, x, fase), T, x)
    return r['h'], r['s']


def estado(T, P, w):
    """Estado completo a (T,P) con composicion global MASICA w."""
    x = ng.w2m(w)
    d = dict(T=T, P=P, w=w)
    eq = flash_TP(T, P)
    if eq is None:
        Ta, Tw = Tsat_puros(P)
        if Ta + 0.3 < T < Tw - 0.3:
            # El flash a (T,P) no converge: ocurre solo en la vecindad de los
            # componentes puros, donde la ventana bifasica en composicion se
            # cierra. Se decide la fase con las funciones propias del motor,
            # que si son estables ahi porque estan parametrizadas por x.
            Pb = ng.bubbleP(T, x)[0]
            Pd = ng.dewP(T, x)[0]
            if P >= Pb:
                fase = 'liquido'
            elif P <= Pd:
                fase = 'vapor'
            else:
                raise RuntimeError(f"equilibrio no resoluble en T={T:.4f} K, "
                                   f"P={P:.5g} MPa, w={w:.4f}")
        else:
            fase = 'liquido' if T <= Ta else 'vapor'
            if Ta < T <= Ta + 0.3:
                # Sobre Ta la regla de componente puro llama vapor a mezclas liquidas.
                try:
                    if P >= ng.bubbleP(T, x)[0]:
                        fase = 'liquido'
                except Exception:
                    pass
        h, s = _mono(T, P, x, 'l' if fase == 'liquido' else 'v')
        d.update(h=h, s=s, q=0.0 if fase == 'liquido' else 1.0, fase=fase)
        return d
    xL, xV = eq
    if x <= xL:
        h, s = _mono(T, P, x, 'l')
        d.update(h=h, s=s, q=0.0, fase='liquido')
    elif x >= xV:
        h, s = _mono(T, P, x, 'v')
        d.update(h=h, s=s, q=1.0, fase='vapor')
    else:
        wL, wV = ng.m2w(xL), ng.m2w(xV)
        q = (w - wL) / (wV - wL)                 # titulo MASICO (regla de la palanca)
        hL, sL = _mono(T, P, xL, 'l')
        hV, sV = _mono(T, P, xV, 'v')
        d.update(h=q * hV + (1 - q) * hL, s=q * sV + (1 - q) * sL, q=q,
                 fase='bifasico', wL=wL, wV=wV, hL=hL, hV=hV, sL=sL, sV=sV)
    return d


_env = {}


def T_sat(P, w, kind):
    """T de liquido saturado ('L') o de vapor saturado ('V') a (P, w masica).
    Solo se usa para informar y para el criterio F1; se cachea."""
    k = (round(P, 10), round(w, 8), kind)
    if k not in _env:
        x = ng.w2m(w)
        _env[k] = (ng.bubbleT(P, x)[0] if kind == 'L' else ng.dewT(P, x)[0])
    return _env[k]


_ultimaT = {}     # (P, w, propiedad) -> ultima T resuelta, para acotar el brentq


def _newton_salva_T(f, lo, hi, flo, fhi, x0, f0=None):
    """Newton-Raphson salvaguardado para invertir h o s -> T (estado_de).

    La derivada se estima en diferencias hacia adelante con paso relativo:
    una sola llamada extra a f() (= estado) por iteracion. Red de seguridad:
    paso fuera del bracket, salto desproporcionado, residuo que no mejora o
    valor no finito (p. ej. cruce de cambio de fase) derivan a biseccion por
    signo; agotado el limite de Newton+bisecciones, brentq es el ultimo
    recurso. Devuelve la raiz con |paso| <= TOL_T (misma tolerancia que el
    brentq que reemplaza)."""
    for _ in range(50):
        if f0 is None:
            f0 = f(x0)
        if not np.isfinite(f0):
            return brentq(f, lo, hi, xtol=TOL_T)
        d = 1e-5 * max(abs(x0), 1.0)
        f1 = f(x0 + d)
        if np.isfinite(f1):
            den = f1 - f0
            if den != 0.0 and np.isfinite(den):
                x1 = x0 - f0 * d / den
                if lo < x1 < hi and abs(x1 - x0) <= 0.5 * (hi - lo):
                    if abs(x1 - x0) <= TOL_T:
                        return x1
                    g1 = f(x1)
                    if np.isfinite(g1) and abs(g1) < abs(f0):
                        x0, f0 = x1, g1
                        continue
        # Biseccion por signo (red de seguridad).
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if not np.isfinite(fm):
            return brentq(f, lo, hi, xtol=TOL_T)
        if fm == 0.0:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if abs(flo) <= abs(fhi):
            x0, f0 = lo, flo
        else:
            x0, f0 = hi, fhi
    return brentq(f, lo, hi, xtol=TOL_T)


def estado_de(P, w, h=None, s=None):
    """Inversion: estado a (P, w) con h o s prescrita. h y s crecen con T.

    El intervalo de arranque se toma alrededor de la ultima T resuelta para esa
    misma (P, w, propiedad) y se ensancha hasta encontrar cambio de signo. En un
    lazo iterativo las llamadas sucesivas caen casi encima, de modo que brentq
    trabaja sobre un intervalo de pocos K en lugar de sobre toda la campana.
    """
    nombre = 'h' if h is not None else 's'
    val = h if h is not None else s
    k = (round(P, 10), round(w, 8), nombre)
    Ta, Tw = Tsat_puros(P)
    T0 = _ultimaT.get(k)
    f = lambda T: estado(T, P, w)[nombre] - val
    if T0 is not None:
        lo, hi = T0 - 8.0, T0 + 8.0
    else:
        lo, hi = max(235.0, Ta - 60.0), Tw + 40.0
    flo, fhi = f(lo), f(hi)
    paso, n = 16.0, 0
    while flo > 0 and n < 30:
        hi, fhi = lo, flo
        lo = max(233.0, lo - paso)
        flo = f(lo)
        paso *= 1.6
        n += 1
        if lo <= 233.0 and flo > 0:
            break
    paso = 16.0
    while fhi < 0 and n < 60:
        lo, flo = hi, fhi
        hi += paso
        fhi = f(hi)
        paso *= 1.6
        n += 1
        if hi > 900.0:
            break
    if flo > 0 or fhi < 0:
        raise RuntimeError(f"{nombre}={val:.6g} fuera de dominio a "
                           f"P={P:.5g} MPa, w={w:.4f}")
    if T0 is not None and lo <= T0 <= hi:
        x0, f0 = T0, None
    else:
        x0, f0 = (lo, flo) if abs(flo) <= abs(fhi) else (hi, fhi)
    T = _newton_salva_T(f, lo, hi, flo, fhi, x0, f0)
    _ultimaT[k] = T
    return estado(T, P, w)


# ---------------------------------------------------------------------------
# CAPA DE RESOLUCION  —  UNICA E IDENTICA PARA LOS DOS METODOS
#
# Los metodos A y B difieren SOLO en las tres funciones de cierre (HRVG,
# regenerador, condensador). Todo lo demas —incluida la rama degenerada y la
# estrategia de solucion— es un solo camino de codigo.
#
# Estructura del sistema:
#   * corriente de corte EXTERIOR: h1  ->  F(h1) = G(h1) - h1 = 0, por Brent.
#     La ausencia de cambio de signo en el intervalo ES la definicion operativa
#     de NO_CONVERGE.
#   * lazo INTERIOR frio (h10 -> h6 -> h7 -> h8 -> h9 -> h10): aparece cuando el
#     condensador cierra por efectividad, porque h9 pasa a depender de h8. Su
#     ganancia es (1-eps_cond)*eps_reg*(m_p/m_b) ~ 0.14: se resuelve por
#     sustitucion sucesiva, que converge monotona en ~10 pasadas.
#     Con el condensador cerrado por approach (T9 dato) la ganancia es cero y
#     el mismo lazo termina en una pasada. No hay bifurcacion de codigo.
#   * en el metodo A, T2 es dato, luego h2 no depende de h1 y G(h1) es constante.
#     El resolvedor lo detecta evaluando G en los dos extremos del intervalo y
#     resuelve sin iterar. Es una propiedad de la fisica, no un atajo.
# ---------------------------------------------------------------------------
MAX_FRIO = 300
# Correccion (2026-09-10, diagnostico de la no-convergencia en la frontera de
# fase): TOL_FRIO estaba en 1e-8, POR DEBAJO del piso de ruido de TOL_T
# (inversion h->T, 1e-6 K sobre cp~4 kJ/kg-K son ~4e-6 kJ/kg) -- el mismo
# patron de las correcciones C3/C23/C23bis, sin corregir aca. El chequeo
# primario (d < TOL_FRIO) nunca disparaba, y CADA llamada a lazo_frio pagaba
# el fallback de 15 iteraciones de estancamiento como piso, no como excepcion.
# Es barato cuando estado_de() esta bien cacheado (fraccion de ms, caso
# degenerado con wp=w_b), pero con una composicion wp nueva (rama Kalina real,
# nunca vista antes por el proceso) cada estado_de() cuesta ~2-4 s, y 15
# iteraciones de mas por CADA evaluacion del lazo exterior son minutos.
# TOL_FRIO ahora esta por ENCIMA del piso de ruido: el camino normal sale por
# el chequeo primario, sin esperar al fallback. TOL_FRIO_PISO se mantiene como
# red de seguridad mas laxa para casos genuinamente lentos en converger, y el
# contador de estancamiento se acorto de 15 a 6 (ya no hace falta tanta
# insistencia si el umbral primario es realista).
TOL_FRIO = 2e-5
TOL_FRIO_PISO = 1e-4
N_ESTANCADO = 6
TOL_H1 = 1e-4             # kJ/kg, ancho final del intervalo de Brent (N1). El
                          # piso de ruido lo pone la inversion h->T (TOL_T):
                          # 1e-6 K sobre cp~4 kJ/kg-K son ~4e-6 kJ/kg. Pedir
                          # 1e-8 (como se hizo antes) repite el error de C3:
                          # doce cifras que el motor no tiene, a cambio de
                          # ~15 iteraciones extra por corrida (correccion C23bis).


class NoConverge(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# LIMITADOR DE PINZAMIENTO POR CAPACIDAD CALORIFICA FINITA DE LA FUENTE (C_g)
#
# Comun a los metodos A y B, se aplica DESPUES del cierre del HRVG (approach o
# efectividad, segun el metodo). Sustituye al factor min(1, phi) del primer
# intento (planteamiento eps-NTU), que es correcto pero inoperante en una
# mezcla zeotropica: lo que muerde no es la capacidad TOTAL del gas sino la
# FORMA del perfil (ver incertidumbre.py, docstring del modulo).
#
# CONVENCION DE FLUJO: CONTRACORRIENTE. El gas entra caliente por el extremo
# T2 (el mas caliente del fluido de trabajo, junto al cierre approach/eps) y
# sale por el extremo T1 a T_gas,out = Tf - Qi/C_g. Es la disposicion estandar
# de cualquier caldera de recuperacion (economizador/evaporador/sobrecalenta-
# dor: el gas mas caliente se encuentra con el fluido mas caliente) y es la
# unica consistente con el cierre por approach ya existente, dT_app = Tf - T2
# (que ya supone que T2, no T1, es el que se acerca a Tf).
#   T_gas(h) = Tf - (h2 - h) / C_g
# ---------------------------------------------------------------------------
def perfil_HRVG(P, w, T1, h1, T2, h2, N=40):
    """Perfil (T,h) del fluido de trabajo en el HRVG entre los estados 1 y 2.

    Se recorre en TEMPERATURA: h(T) es evaluacion directa, T(h) exigiria una
    inversion por nodo. Los extremos se fuerzan a los valores REALES del
    cierre (T1,h1) y (T2,h2) en vez de dejar que caigan donde toque la
    rejilla uniforme: sin esto el minimo del extremo caliente -justo donde
    vive el approach declarado- se pierde por un artefacto de rejilla. El
    perfil se autovalida asi contra el cierre que lo genero.
    """
    if T2 <= T1:
        return np.array([T1, T2]), np.array([h1, h2])
    Ts = np.linspace(T1, T2, N)
    hs = np.array([estado(float(T), P, w)['h'] for T in Ts])
    Ts[0], hs[0] = T1, h1
    Ts[-1], hs[-1] = T2, h2
    return Ts, hs


def pinzamiento(P, w, T1, h1, T2, h2, Tf, C_g, dT_pp):
    """Recorta (h2, T2) si el margen dT_pp no se cumple en algun punto del
    HRVG con esa capacidad calorifica de fuente C_g (contracorriente).

    Estrategia: tabular el perfil UNA vez sobre el rango declarado [T1,T2] (N
    evaluaciones de estado()) y resolver el recorte por interpolacion sobre esa
    misma tabla -sin llamadas nuevas a estado()-, para no reintroducir el coste
    de un lazo interno de Brent por cada llamada. Devuelve dict(h2, T2, DT_min,
    z_min, recortado).
    """
    if T2 <= T1 or not np.isfinite(C_g):
        return dict(h2=h2, T2=T2, DT_min=np.inf, z_min=0.0, recortado=False)
    Ts, hs = perfil_HRVG(P, w, T1, h1, T2, h2)
    Tg = Tf - (h2 - hs) / C_g
    DT = Tg - Ts
    j = int(np.argmin(DT))
    den = hs[-1] - hs[0]
    z = (hs - hs[0]) / den if den > 0 else np.zeros_like(hs)
    if DT[j] >= dT_pp - 1e-9:
        return dict(h2=h2, T2=T2, DT_min=float(DT[j]), z_min=float(z[j]),
                    recortado=False)
    # bound[i]: valor de h2 que deja el nodo i exactamente al margen dT_pp,
    # tratando h2 como variable libre e i como nodo interior fijo:
    #   Tf - (h2-h_i)/C_g - T_i >= dT_pp  <=>  h2 <= h_i + C_g*(Tf-T_i-dT_pp)
    # El h2 factible mas grande es el menor de estos bounds sobre TODOS los
    # nodos que quedarian dentro del equipo para ese h2 (running-min desde el
    # extremo frio, que es donde arranca cualquier cierre valido).
    bound = hs + C_g * (Tf - Ts - dT_pp)
    run_min = np.minimum.accumulate(bound)
    malo = np.where(hs > run_min + 1e-12)[0]
    k = int(malo[0]) if len(malo) else len(hs)
    if k == 0:
        return dict(h2=hs[0], T2=Ts[0], DT_min=dT_pp, z_min=0.0, recortado=True)
    if k >= len(hs):
        k = len(hs) - 1
    lim = run_min[k - 1]
    f = float(np.clip((lim - hs[k - 1]) / (hs[k] - hs[k - 1]), 0.0, 1.0))
    h2r = hs[k - 1] + f * (hs[k] - hs[k - 1])
    T2r = Ts[k - 1] + f * (Ts[k] - Ts[k - 1])
    zr = (h2r - hs[0]) / den if den > 0 else 0.0
    return dict(h2=float(h2r), T2=float(T2r), DT_min=dT_pp, z_min=float(zr),
                recortado=True)


def C_g_minimo_HRVG(Ts, hs, Tf, dT_pp=0.0):
    """Menor C_g (contracorriente) que evita DT < dT_pp en todo el HRVG con el
    (T1,h1)-(T2,h2) YA resueltos. Informativo: se reporta en S9 cuando C_g no
    se declara, para dar la capacidad minima de fuente que el punto exige."""
    h2 = hs[-1]
    den = Tf - Ts - dT_pp
    # En el ultimo nodo (h_i=h2) el numerador es 0: si el propio cierre ya
    # incumple dT_pp en el extremo caliente (den<0 ahi), NINGUN C_g finito lo
    # arregla, y el np.where debe propagar ese inf, no ocultarlo.
    req = np.where(den > 0, (h2 - hs) / np.maximum(den, 1e-300), np.inf)
    req = np.where((h2 - hs) <= 0, 0.0, req)
    return float(np.max(req)) if len(req) else 0.0


def _cg_min_hrvg(est, ind, par):
    """C_g,min del HRVG (kW/K por kg/s de m_b) bajo demanda, sobre la solucion
    ya resuelta en est/ind y los parametros par. Lo consume solo el TEXTO de
    supuestos()/criterios(), por eso no esta en el camino caliente de
    resolver(). Devuelve None si el perfil no se puede evaluar."""
    try:
        Ts, hs = perfil_HRVG(par['P_alta'], par['w_b'], est[1]['T'], est[1]['h'],
                             est[2]['T'], est[2]['h'])
        return C_g_minimo_HRVG(Ts, hs, par['T_f'], par.get('dT_pp_gas', 0.0))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# CONDENSADOR CON AGUA DE CAPACIDAD CALORIFICA FINITA (extension C_cold)
#
# Espejo del limitador del HRVG, del lado frio. El agua entra a T_agua_in por el
# extremo del estado 9 y sale por el del estado 8 (contracorriente):
#   T_w(h) = T_agua_in + (h - h9)/C_cold,  h en [h9, h8]
# La mezcla condensa con deslizamiento, asi que el acercamiento minimo puede
# quedar DENTRO del equipo y no en el extremo frio donde se declara el cierre.
# C_cold no declarado reproduce EXACTAMENTE el sumidero infinito anterior.
# ---------------------------------------------------------------------------
DT_PERFIL_COND = 2.0      # paso en temperatura del perfil del condensador [K]
SPAN_PERFIL_COND = 300.0  # tramo maximo explorado por encima de T9 [K]


def perfil_condensador(P, w, T9, h9, h8, dT=DT_PERFIL_COND, span=SPAN_PERFIL_COND):
    """Perfil (T,h) de la solucion basica en el condensador, de 9 a 8.

    Correccion C31: NO invierte T8 = T(h8). Esa inversion, repetida en cada
    pasada del lazo frio, caia en la zona donde el flash del motor es fragil
    (cerca del rocio a P_baja) y el ciclo no convergia. Aqui se sube en
    temperatura desde T9 con paso dT, con evaluaciones directas h(T); un nodo
    donde estado() falla se SALTA y se cuenta. T8 se obtiene interpolando entre
    el ultimo nodo valido y el primero que supera h8. Devuelve (Ts, hs,
    nodos_saltados); el ultimo nodo es exactamente (T8, h8).
    """
    Ts, hs = [T9], [h9]
    if h8 <= h9:
        return np.array(Ts), np.array(hs), 0
    saltados = 0
    T = T9
    while T < T9 + span:
        T += dT
        try:
            h = estado(T, P, w)['h']
        except RuntimeError:
            saltados += 1
            continue
        if h >= h8:
            Tp, hp = Ts[-1], hs[-1]
            Ts.append(Tp + (T - Tp) * (h8 - hp) / (h - hp))
            hs.append(h8)
            return np.array(Ts), np.array(hs), saltados
        if h > hs[-1]:
            Ts.append(T)
            hs.append(h)
    raise RuntimeError(f"perfil del condensador: h8 no alcanzado en {span:.0f} K "
                       f"sobre T9={T9:.2f} K (P={P:.5g} MPa, w={w:.4f})")


def pinzamiento_condensador(Ts, hs, T_w_in, C_cold, dT_pp):
    """Sube (h9, T9) -- menos calor retirado -- si el margen dT_pp no se cumple
    en algun punto del condensador con esa capacidad de agua C_cold. Mismo
    esquema que pinzamiento(): sobre la tabla (Ts, hs) de perfil_condensador(),
    recorte por interpolacion. Devuelve dict(h9, T9, DT_min, z_min, recortado,
    sin_solucion).
    """
    T9, h9 = float(Ts[0]), float(hs[0])
    if len(hs) < 2 or not np.isfinite(C_cold):
        return dict(h9=h9, T9=T9, DT_min=np.inf, z_min=0.0, recortado=False,
                    sin_solucion=False)
    DT = Ts - (T_w_in + (hs - h9) / C_cold)
    j = int(np.argmin(DT))
    den = hs[-1] - hs[0]
    z = (hs - hs[0]) / den if den > 0 else np.zeros_like(hs)
    if DT[j] >= dT_pp - 1e-9:
        return dict(h9=h9, T9=T9, DT_min=float(DT[j]), z_min=float(z[j]),
                    recortado=False, sin_solucion=False)
    # bound[i]: menor h9 que deja el nodo i exactamente al margen dT_pp:
    #   T_i - T_w_in - (h_i-h9)/C >= dT_pp  <=>  h9 >= h_i - C*(T_i-T_w_in-dT_pp)
    # El h9 factible mas chico debe superar los bounds de TODOS los nodos que
    # quedan dentro del equipo (los de h >= h9): maximo acumulado desde el
    # extremo caliente. g = h - run_max es creciente, se busca su cero.
    bound = hs - C_cold * (Ts - T_w_in - dT_pp)
    run_max = np.maximum.accumulate(bound[::-1])[::-1]
    g = hs - run_max
    buenos = np.where(g >= -1e-12)[0]
    if len(buenos) == 0:
        return dict(h9=float(hs[-1]), T9=float(Ts[-1]), DT_min=-np.inf, z_min=1.0,
                    recortado=True, sin_solucion=True)
    k = max(int(buenos[0]), 1)
    f = float(np.clip(-g[k - 1] / (g[k] - g[k - 1]), 0.0, 1.0))
    h9r = hs[k - 1] + f * (hs[k] - hs[k - 1])
    T9r = Ts[k - 1] + f * (Ts[k] - Ts[k - 1])
    zr = (h9r - hs[0]) / den if den > 0 else 0.0
    return dict(h9=float(h9r), T9=float(T9r), DT_min=dT_pp, z_min=float(zr),
                recortado=True, sin_solucion=False)


def C_cold_minimo(Ts, hs, T_w_in, dT_pp=0.0):
    """Menor C_cold (contracorriente) que evita DT < dT_pp en todo el
    condensador con el (T9,h9)-(T8,h8) YA resueltos. Informativo: es la
    condicion que el resultado exige al agua de enfriamiento cuando C_cold no
    se declara. inf si el propio extremo frio ya incumple el margen."""
    if len(hs) == 0:
        return 0.0
    h9 = hs[0]
    den = Ts - T_w_in - dT_pp
    if den[0] < 0:
        return float('inf')
    req = np.where(den > 0, (hs - h9) / np.maximum(den, 1e-300), np.inf)
    req = np.where((hs - h9) <= 0, 0.0, req)
    return float(np.max(req))


def resolver(par, devolver_traza=False):
    """Resuelve los 10 estados del KCS-11.

    par:
      P_alta, P_baja [MPa]; w_b [masica]; T_f, T_amb [K]; eta_t, eta_p
      cierre_HRVG = {'dT_app': K}  (A)  |  {'eps': -}  (B)
      cierre_reg  = {'dT_pp': K}   (A)  |  {'eps': -}  (B)
      cierre_cond = {'dT_app': K}  (A)  |  {'eps': -}  (B)

    Opcionales (si faltan se asume la hipotesis indicada y supuestos() declara
    la incertidumbre que eso deja):
      C_g [kW/K por kg/s de m_b]    fuente           -> reservorio isotermo (H10)
      dT_pp_gas [K]                 margen HRVG      -> 0
      T_min_gas [K]                 rocio acido (O6) -> no evaluado
      C_cold [kW/K por kg/s de m_b] agua de enfriam. -> sumidero infinito
      T_agua_in [K]                 entrada del agua -> T_amb
      dT_pp_cond [K]                margen condens.  -> 0
      T_agua_out_max [K]            salida agua (O7) -> no evaluado
      cp_agua [kJ/kg-K]             solo para reportar caudales de agua -> 4.18
    """
    Ph, Pl = par['P_alta'], par['P_baja']
    wb, Tf, Ta = par['w_b'], par['T_f'], par['T_amb']
    et, ep = par['eta_t'], par['eta_p']
    cH, cR, cC = par['cierre_HRVG'], par['cierre_reg'], par['cierre_cond']
    # Temperatura del sumidero que ve el condensador. Sin T_agua_in es T_amb,
    # que es exactamente el comportamiento anterior. T_amb sigue siendo el
    # estado muerto (Carnot, eta_recurso).
    Tk = par.get('T_agua_in', Ta)

    h2_ideal = estado(Tf, Ph, wb)['h']          # equilibrio termico con la fuente
    h9_ideal = estado(Tk, Pl, wb)['h']          # equilibrio termico con el sumidero

    # ---- cierres --------------------------------------------------------
    # C_g [kW/K por kg/s de m_b]: capacidad calorifica de la fuente, declarada
    # de forma ABSOLUTA (no como razon phi=C_g/C_wf: C_wf no se conoce antes de
    # resolver, y el planteamiento eps-NTU basado en phi resulto INOPERANTE en
    # una mezcla zeotropica -- ver docstring de incertidumbre.py). C_g=None
    # (o ausente) reproduce EXACTAMENTE la hipotesis H10 (fuente como
    # reservorio) y deja el limitador de pinzamiento inactivo.
    C_g = par.get('C_g')
    dT_pp_gas = par.get('dT_pp_gas', 0.0)
    pinch_info = dict(DT_min=np.inf, z_min=None, recortado=False, Cg_min=None)
    # Sonda de incertidumbre de composicion: desplaza w_r y w_p dentro de la banda
    # declarada por G4-01 (+-0.01). NO es un estado consistente: se reporta el
    # residuo del balance de energia del separador como medida del desvio.
    dwr, dwp = par.get('perturb_w', (0.0, 0.0))
    # C_cold: mismo convenio que C_g. None (o ausente) = sumidero infinito,
    # limitador del condensador inactivo, comportamiento anterior intacto.
    C_cold = par.get('C_cold')
    dT_pp_cond = par.get('dT_pp_cond', 0.0)
    cold_activo = C_cold is not None and np.isfinite(C_cold)
    cond_info = dict(DT_min=np.inf, z_min=None, recortado=False, sin_solucion=False,
                     nodos_saltados=0)

    _debug = par.get('_debug')
    _dbg_n = [0]

    def cierre_hrvg(h1, T1):
        if 'dT_app' in cH:
            e2 = estado(Tf - cH['dT_app'], Ph, wb)
        else:
            e2 = estado_de(Ph, wb, h=h1 + cH['eps'] * (h2_ideal - h1))
        if C_g is not None and np.isfinite(C_g):
            r = pinzamiento(Ph, wb, T1, h1, e2['T'], e2['h'], Tf, C_g, dT_pp_gas)
            pinch_info.update(DT_min=r['DT_min'], z_min=r['z_min'],
                               recortado=r['recortado'])
            if r['recortado']:
                e2 = estado(r['T2'], Ph, wb)
        if _debug:
            _dbg_n[0] += 1
            print(f"[{_dbg_n[0]:4d}] {time.time():.2f}  T1={T1:.4f}  "
                  f"T2_decl={e2['T']:.4f} fase={e2['fase']:8s} q2={e2['q']:.4f} "
                  f"recorte={pinch_info['recortado']}", flush=True)
        return e2

    def cierre_reg(T10, wp, h5):
        h6_ideal = estado(T10, Ph, wp)['h']
        if 'dT_pp' in cR:
            return estado(T10 + cR['dT_pp'], Ph, wp), h6_ideal
        return estado_de(Ph, wp, h=h5 - cR['eps'] * (h5 - h6_ideal)), h6_ideal

    def cierre_cond(h8):
        if 'dT_app' in cC:
            e9 = estado(Tk + cC['dT_app'], Pl, wb)
        else:
            e9 = estado_de(Pl, wb, h=h8 - cC['eps'] * (h8 - h9_ideal))
        if cold_activo:
            Ts_c, hs_c, n_salt = perfil_condensador(Pl, wb, e9['T'], e9['h'], h8)
            r = pinzamiento_condensador(Ts_c, hs_c, Tk, C_cold, dT_pp_cond)
            cond_info.update(DT_min=r['DT_min'], z_min=r['z_min'],
                             recortado=r['recortado'],
                             sin_solucion=r['sin_solucion'],
                             nodos_saltados=n_salt)
            if r['recortado']:
                e9 = estado(r['T9'], Pl, wb)
        return e9

    # ---- lazo interior frio --------------------------------------------
    def lazo_frio(h4, mr, mp, wp, h5):
        h10 = h9_ideal + 5.0
        e6 = e7 = e9 = e10 = None
        d_min, estancado = np.inf, 0
        for it in range(MAX_FRIO):
            if _debug:
                _t0 = time.time()
                e10 = estado_de(Ph, wb, h=h10)
                _t1 = time.time()
                e6, h6_ideal = cierre_reg(e10['T'], wp, h5)
                _t2 = time.time()
                e7 = estado_de(Pl, wp, h=e6['h']) if mp > 0 else dict(
                    T=e6['T'], P=Pl, w=wp, h=e6['h'], s=e6['s'], q=0.0, fase='liquido')
                _t3 = time.time()
                h8 = mr * h4 + mp * e7['h']
                e9 = cierre_cond(h8)
                _t4 = time.time()
                e10s = estado_de(Ph, wb, s=e9['s'])
                _t5 = time.time()
                print(f"    frio[{it:3d}] e10={_t1-_t0:6.2f}s e6={_t2-_t1:6.2f}s "
                      f"e7={_t3-_t2:6.2f}s e9={_t4-_t3:6.2f}s e10s={_t5-_t4:6.2f}s "
                      f"wp={wp:.5f} h6={e6['h']:.2f}", flush=True)
            else:
                e10 = estado_de(Ph, wb, h=h10)
                e6, h6_ideal = cierre_reg(e10['T'], wp, h5)
                e7 = estado_de(Pl, wp, h=e6['h']) if mp > 0 else dict(
                    T=e6['T'], P=Pl, w=wp, h=e6['h'], s=e6['s'], q=0.0, fase='liquido')
                h8 = mr * h4 + mp * e7['h']
                e9 = cierre_cond(h8)
                e10s = estado_de(Ph, wb, s=e9['s'])
            h10n = e9['h'] + (e10s['h'] - e9['h']) / ep
            d = abs(h10n - h10)
            h10 = h10n
            # Estancamiento en el piso de ruido: la sustitucion sucesiva ya no
            # mejora porque la inversion h->T no da mas cifras. Se acepta.
            if d < d_min - 1e-12:
                d_min, estancado = d, 0
            else:
                estancado += 1
            if d < TOL_FRIO or (estancado >= N_ESTANCADO and d < TOL_FRIO_PISO):
                e10 = estado_de(Ph, wb, h=h10)
                e6, h6_ideal = cierre_reg(e10['T'], wp, h5)
                e7 = estado_de(Pl, wp, h=e6['h']) if mp > 0 else e7
                h8 = mr * h4 + mp * e7['h']
                e9 = cierre_cond(h8)
                return dict(e6=e6, e7=e7, e9=e9, e10=e10, h8=h8, h6_ideal=h6_ideal,
                            e10s=e10s, iter=it + 1)
        raise NoConverge(f"lazo frio no converge (residuo {d:.3e} kJ/kg)")

    # ---- rama alta: de h1 a h4, y de ahi al lazo frio -------------------
    def tramo(h1):
        T1 = estado_de(Ph, wb, h=h1)['T']
        e2 = cierre_hrvg(h1, T1)
        if e2['fase'] == 'bifasico':
            wr, wp = e2['wV'] + dwr, e2['wL'] + dwp
            mr = (wb - wp) / (wr - wp)
            mp = 1.0 - mr
            if dwr or dwp:      # sonda de incertidumbre: estados recalculados
                e3, e5 = estado(e2['T'], Ph, wr), estado(e2['T'], Ph, wp)
            else:
                e3 = dict(T=e2['T'], P=Ph, w=wr, h=e2['hV'], s=e2['sV'],
                          q=1.0, fase='vapor')
                e5 = dict(T=e2['T'], P=Ph, w=wp, h=e2['hL'], s=e2['sL'],
                          q=0.0, fase='liquido')
            degenerado = False
        else:
            # RAMA DEGENERADA: el separador no separa. m_p = 0, m_r = m_b,
            # x_r = x_b. Mismo camino de codigo; solo cambian tres asignaciones.
            wr, wp = wb, wb
            mr, mp = 1.0, 0.0
            e3 = dict(e2)
            e5 = dict(e2)
            degenerado = True
        e4s = estado_de(Pl, wr, s=e3['s'])
        h4 = e3['h'] - et * (e3['h'] - e4s['h'])
        e4 = estado_de(Pl, wr, h=h4)
        frio = lazo_frio(h4, mr, mp, wp, e5['h'])
        h1_nuevo = frio['e10']['h'] + mp * (e5['h'] - frio['e6']['h'])
        return h1_nuevo, dict(e2=e2, e3=e3, e4=e4, e4s=e4s, e5=e5, wr=wr, wp=wp,
                              mr=mr, mp=mp, degenerado=degenerado, **frio)

    # ---- lazo exterior: F(h1) = G(h1) - h1 = 0 --------------------------
    h1_lo_full = estado(Tk, Ph, wb)['h']       # cota inferior: sin regeneracion
    h1_hi_full = h2_ideal                      # cota superior: equilibrio con la fuente
    # Arranque tibio: si el llamador pasa par['h1_semilla'] (tipicamente el h1
    # resuelto en el punto vecino de un barrido), se acota el intervalo de
    # Brent alrededor de esa semilla en vez de partir del rango fisico
    # completo, y se ensancha solo si no hay cambio de signo ahi. Mismo
    # principio que la memoria de intervalo de estado_de(). Sin semilla
    # (par no la trae, caso por defecto: prueba4/5 y cualquier corrida
    # aislada) el comportamiento es EXACTAMENTE el de antes.
    semilla = par.get('h1_semilla')
    if semilla is None:
        h1_lo, h1_hi = h1_lo_full, h1_hi_full
        G_lo, _ = tramo(h1_lo)
        G_hi, _ = tramo(h1_hi)
    else:
        semilla = min(max(semilla, h1_lo_full), h1_hi_full)
        ancho_s = par.get('h1_semilla_ancho', 2.0)
        n = 0
        while True:
            h1_lo = max(h1_lo_full, semilla - ancho_s)
            h1_hi = min(h1_hi_full, semilla + ancho_s)
            G_lo, _ = tramo(h1_lo)
            G_hi, _ = tramo(h1_hi)
            bracket_ok = (G_lo - h1_lo) * (G_hi - h1_hi) <= 0
            en_tope = h1_lo <= h1_lo_full + 1e-9 and h1_hi >= h1_hi_full - 1e-9
            if bracket_ok or en_tope or n >= 15:
                break
            ancho_s *= 3.0
            n += 1
    if abs(G_hi - G_lo) <= 1e-9:               # metodo A: G constante
        h1 = G_lo
        ancho = 0.0
    else:
        F = lambda h: tramo(h)[0] - h
        if (G_lo - h1_lo) * (G_hi - h1_hi) > 0:
            raise NoConverge("F(h1) no cambia de signo en [h1_lo, h1_hi]")
        h1 = brentq(F, h1_lo, h1_hi, xtol=TOL_H1)
        ancho = TOL_H1
    _, r = tramo(h1)

    e1 = estado_de(Ph, wb, h=h1)
    e2, e3, e4, e5, e6, e7, e9, e10 = (r['e2'], r['e3'], r['e4'], r['e5'],
                                       r['e6'], r['e7'], r['e9'], r['e10'])
    e8 = estado_de(Pl, wb, h=r['h8'])
    mr, mp, wr, wp = r['mr'], r['mp'], r['wr'], r['wp']

    Wp = e10['h'] - e9['h']
    Wt = mr * (e3['h'] - e4['h'])
    Wn = Wt - Wp
    Qi = e2['h'] - h1
    Qo = r['h8'] - e9['h']
    eta_c = 1 - Ta / Tf

    est = {1: e1, 2: e2, 3: e3, 4: e4, 5: e5, 6: e6, 7: e7, 8: e8, 9: e9, 10: e10}
    ind = dict(Qi=Qi, Qout=Qo, Wt=Wt, Wp=Wp, Wnet=Wn, eta=Wn / Qi, h1=h1,
               eta_carnot=eta_c, eta_II=(Wn / Qi) / eta_c, RBW=Wp / Wt,
               m_r=mr, m_p=mp, omega=mp, w_r=wr, w_p=wp,
               q2=e2['q'], q4=e4['q'], q7=e7['q'], q8=e8['q'], q9=e9['q'],
               eps_HRVG=(e2['h'] - h1) / (h2_ideal - h1),
               dT_app=Tf - e2['T'],
               eps_reg=((e5['h'] - e6['h']) / (e5['h'] - r['h6_ideal'])
                        if mp > 0 else float('nan')),
               dT_pp=e6['T'] - e10['T'],
               eps_cond=(r['h8'] - e9['h']) / (r['h8'] - h9_ideal),
               dT_app_cond=e9['T'] - Tk,
               degenerado=r['degenerado'], iter_frio=r['iter'], ancho_h1=ancho,
               T1=e1['T'], T2=e2['T'], T9=e9['T'], T10=e10['T'])

    # --- capacidad calorifica de la fuente y sus consecuencias --------------
    C_wf = (h2_ideal - h1) / (Tf - e1['T'])        # kW/K por kg/s de fluido basico
    ind['C_wf'] = C_wf
    ind['res_separador'] = e2['h'] - (mr * e3['h'] + mp * e5['h'])
    ind['DT_min_HRVG'] = pinch_info['DT_min']
    ind['z_min_HRVG'] = pinch_info['z_min']
    ind['pinch_recortado'] = pinch_info['recortado']
    if C_g is not None and np.isfinite(C_g):
        ind['C_g'] = C_g
        ind['phi'] = C_g / C_wf
        ind['T_gas_out'] = Tf - Qi / C_g
        ind['eta_recurso'] = Wn / (C_g * (Tf - Ta))
        ind['Cg_min_HRVG'] = None
    else:
        ind['C_g'] = np.inf
        ind['phi'] = np.inf
        ind['T_gas_out'] = Tf
        ind['eta_recurso'] = 0.0
        # C_g,min es solo INFORMATIVO (lo leen supuestos()/criterios() para
        # armar texto): se reconstruye bajo demanda con _cg_min_hrvg() y no
        # se paga el perfil en el camino caliente de resolver().
        ind['Cg_min_HRVG'] = None

    # --- agua de enfriamiento y sus consecuencias ---------------------------
    cp_w = par.get('cp_agua', 4.18)
    ind['T_agua_in'] = Tk
    ind['DT_min_cond'] = cond_info['DT_min']
    ind['z_min_cond'] = cond_info['z_min']
    ind['cond_recortado'] = cond_info['recortado']
    ind['cond_sin_solucion'] = cond_info['sin_solucion']
    ind['cond_nodos_saltados'] = cond_info['nodos_saltados']
    if cold_activo:
        ind['C_cold'] = C_cold
        ind['T_agua_out'] = Tk + Qo / C_cold
        ind['m_agua'] = C_cold / cp_w
    else:
        ind['C_cold'] = np.inf
        ind['T_agua_out'] = Tk
        ind['m_agua'] = np.inf
    # C_cold_min NO se calcula aqui: evaluar el perfil llena las memorias de
    # arranque tibio del motor y podria mover en el ruido las corridas que no
    # declaran nada. Se calcula bajo demanda en diagnostico_condensador().
    if devolver_traza:
        return est, ind, r
    return est, ind


# ---------------------------------------------------------------------------
# CRITERIOS DE ADMISIBILIDAD (subconjunto activo en la prueba 4)
# ---------------------------------------------------------------------------
def criterios(est, ind, par):
    e, wb, Ph = est, par['w_b'], par['P_alta']
    deg = ind['degenerado']
    c = []
    add = lambda i, ok, sev, txt: c.append((i, bool(ok), sev, txt))

    r = abs(ind['Qi'] - ind['Qout'] - ind['Wnet'])
    add('N2', r <= 1e-6 * ind['Qi'], 'numerico',
        f"|Qi-Qout-Wnet| = {r:.3e} kJ/kg   (tol = {1e-6*ind['Qi']:.3e})")
    r1 = abs(wb - ind['m_r'] * ind['w_r'] - ind['m_p'] * ind['w_p'])
    add('N3', r1 <= 1e-6 * wb, 'numerico',
        f"separador y absorbedor: residuo {r1:.2e}")
    add('S1', e[4]['s'] >= e[3]['s'], 'critico',
        f"s4-s3 = {e[4]['s']-e[3]['s']:+.6f} kJ/kg-K")
    add('S2', e[10]['s'] >= e[9]['s'], 'critico',
        f"s10-s9 = {e[10]['s']-e[9]['s']:+.6f} kJ/kg-K")
    # Los estados 6 y 7 son de la corriente pobre, que no lleva masa en la rama degenerada.
    if not deg:
        add('S3', e[7]['s'] >= e[6]['s'], 'critico',
            f"s7-s6 = {e[7]['s']-e[6]['s']:+.6f} kJ/kg-K")
    else:
        add('S3', True, 'no aplica', "rama degenerada: la corriente pobre no lleva masa")
    add('S4', e[5]['T'] > e[1]['T'] and e[6]['T'] > e[10]['T'], 'critico',
        f"T5-T1 = {e[5]['T']-e[1]['T']:+.2f} K ; T6-T10 = {e[6]['T']-e[10]['T']:+.2f} K")
    add('S5', e[2]['T'] <= par['T_f'], 'critico',
        f"T2 = {e[2]['T']:.2f} K <= T_f = {par['T_f']:.2f} K")
    # S9 sustituye a la nota de alcance de S5 (correccion C12 del V3): S5 solo
    # ve los extremos del HRVG, no un cruce interno del perfil (posible en una
    # mezcla zeotropica con fuente de capacidad calorifica finita). Con fuente
    # isoterma (C_g no declarado) es EQUIVALENTE a S5 y se marca "no aplica",
    # reportando C_g,min como dato; con C_g finito es estrictamente mas fuerte,
    # y por construccion del limitador de pinzamiento (pinzamiento() recorta
    # h2 hasta cumplirlo) debiera darse siempre por satisfecho -- es tambien
    # una red de consistencia interna del propio resolvedor.
    C_g_ = par.get('C_g')
    dT_pp_gas_ = par.get('dT_pp_gas', 0.0)
    if C_g_ is not None and np.isfinite(C_g_):
        add('S9', ind['DT_min_HRVG'] >= dT_pp_gas_ - 1e-6, 'critico',
            f"DT_min = {ind['DT_min_HRVG']:.3f} K en z = {ind['z_min_HRVG']:.3f} "
            f"(margen declarado {dT_pp_gas_:.1f} K, C_g = {C_g_:.3f} kW/K)"
            + ("  [recorte activo]" if ind['pinch_recortado'] else ""))
    else:
        add('S9', True, 'no aplica',
            f"fuente sin declarar (H10): C_g,min = {_cg_min_hrvg(est, ind, par):.3f} kW/K "
            f"para no cruzar con margen {dT_pp_gas_:.1f} K")
    # O6 -- rocio acido: T_gas,out no puede bajar del punto de rocio acido +
    # margen de diseno. Solo tiene contenido si se declara T_min_gas.
    T_min_gas_ = par.get('T_min_gas')
    if T_min_gas_ is not None:
        add('O6', ind['T_gas_out'] >= T_min_gas_, 'critico',
            f"T_gas,out = {ind['T_gas_out']:.2f} K >= T_min = {T_min_gas_:.2f} K "
            f"(rocio acido + margen)")
    else:
        add('O6', True, 'no aplica', "T_min_gas no declarado")
    # CD -- cierre declarado alcanzable: ADVERTENCIA, no critico. El punto no
    # es imposible (el limitador ya lo resolvio recortando h2); lo que deja de
    # ser cierto es que el cierre DECLARADO (dT_app o eps_HRVG) describe el
    # equipo. Solo tiene contenido si hubo recorte.
    if C_g_ is not None and np.isfinite(C_g_) and ind['pinch_recortado']:
        if 'dT_app' in par['cierre_HRVG']:
            add('CD', False, 'advertencia',
                f"dT_app declarado = {par['cierre_HRVG']['dT_app']:.2f} K; "
                f"realizado = {par['T_f']-e[2]['T']:.2f} K (limitado por C_g)")
        else:
            add('CD', False, 'advertencia',
                f"eps_HRVG declarado = {par['cierre_HRVG']['eps']:.4f}; "
                f"realizado = {ind['eps_HRVG']:.4f} (limitado por C_g)")
    else:
        add('CD', True, 'no aplica',
            "cierre declarado alcanzable, o C_g no declarado")
    # Criterios del agua de enfriamiento. Solo se AGREGAN a la lista si se
    # declara la variable correspondiente: asi la salida de las corridas que no
    # declaran nada (pruebas 4 y 5 incluidas) no cambia. La ausencia se declara
    # en supuestos(), no aqui.
    Tk_ = par.get('T_agua_in', par['T_amb'])
    if 'T_agua_in' in par:
        add('S6', e[9]['T'] >= Tk_, 'critico',
            f"T9-T_agua_in = {e[9]['T']-Tk_:+.2f} K")
    else:
        add('S6', e[9]['T'] >= par['T_amb'], 'critico',
            f"T9-T_amb = {e[9]['T']-par['T_amb']:+.2f} K")
    C_cold_ = par.get('C_cold')
    dT_pp_cond_ = par.get('dT_pp_cond', 0.0)
    if C_cold_ is not None and np.isfinite(C_cold_):
        add('S10', ind['DT_min_cond'] >= dT_pp_cond_ - 1e-6, 'critico',
            ("sin solucion: ningun h9 cumple el margen con ese C_cold"
             if ind['cond_sin_solucion'] else
             f"DT_min = {ind['DT_min_cond']:.3f} K en z = {ind['z_min_cond']:.3f} "
             f"(margen declarado {dT_pp_cond_:.1f} K, C_cold = {C_cold_:.3f} kW/K)")
            + ("  [recorte activo]" if ind['cond_recortado'] else "")
            + (f"  [{ind['cond_nodos_saltados']} nodos sin evaluar]"
               if ind['cond_nodos_saltados'] else ""))
        if ind['cond_recortado']:
            if 'dT_app' in par['cierre_cond']:
                add('CDc', False, 'advertencia',
                    f"dT_app_cond declarado = {par['cierre_cond']['dT_app']:.2f} K; "
                    f"realizado = {ind['dT_app_cond']:.2f} K (limitado por C_cold)")
            else:
                add('CDc', False, 'advertencia',
                    f"eps_cond declarado = {par['cierre_cond']['eps']:.4f}; "
                    f"realizado = {ind['eps_cond']:.4f} (limitado por C_cold)")
    if 'T_agua_out_max' in par:
        if C_cold_ is not None and np.isfinite(C_cold_):
            add('O7', ind['T_agua_out'] <= par['T_agua_out_max'], 'critico',
                f"T_agua,out = {ind['T_agua_out']:.2f} K <= "
                f"T_max = {par['T_agua_out_max']:.2f} K")
        else:
            add('O7', True, 'no aplica',
                "T_agua_out_max declarado pero C_cold no: el agua no se calienta")
    add('S8', ind['eta'] < ind['eta_carnot'], 'critico',
        f"eta = {100*ind['eta']:.2f} % < eta_Carnot = {100*ind['eta_carnot']:.2f} %")
    add('O1', ind['q4'] >= 0.90, 'tecnologico', f"q4 = {ind['q4']:.4f}")
    add('O2', e[9]['fase'] == 'liquido', 'critico',
        f"estado 9 {e[9]['fase']} (q9 = {e[9]['q']:.4f})")
    add('O3', ind['Wnet'] > 0, 'critico', f"W_net = {ind['Wnet']:.3f} kJ/kg")
    add('O5', DQ <= ind['q2'] <= 1 - DQ, 'topologico', f"q2 = {ind['q2']:.4f}")
    # C1 y C1b solo tienen contenido si hay separacion. En la rama degenerada
    # w_r = w_p = w_b por construccion, y evaluarlos como criticos convertiria
    # un DEGENERADO en INVIABLE, que es justo la confusion que el V3 prohibe.
    if not deg:
        add('C1', ind['w_p'] < wb < ind['w_r'], 'critico',
            f"w_p = {ind['w_p']:.4f} < w_b = {wb:.4f} < w_r = {ind['w_r']:.4f}")
        add('C1b', ind['w_r'] - ind['w_p'] > 0.01, 'advertencia',
            f"w_r-w_p = {ind['w_r']-ind['w_p']:.4f}  (banda G4-01: 0.01)")
    else:
        add('C1', True, 'no aplica', "rama degenerada: no hay separacion")
        add('C1b', True, 'no aplica', "rama degenerada: no hay separacion")
    add('C3', 0 <= ind['m_p'] <= 1 and 0 <= ind['m_r'] <= 1, 'critico',
        f"m_r = {ind['m_r']:.4f} ; m_p = {ind['m_p']:.4f}")
    if 'dT_pp' in par['cierre_reg']:
        TL = T_sat(Ph, wb, 'L')
        add('F1', e[1]['T'] < TL, 'advertencia',
            f"T1 = {e[1]['T']:.2f} K vs T_sat,L(P_alta,w_b) = {TL:.2f} K")
        add('PP', (e[5]['T'] - e[1]['T']) >= par['cierre_reg']['dT_pp'] - 1e-9,
            'advertencia',
            f"DT extremo caliente = {e[5]['T']-e[1]['T']:.2f} K "
            f"(pinch declarado {par['cierre_reg']['dT_pp']:.1f} K en el extremo frio)")
    return c


def clasificar(crit):
    """Etiqueta de salida y lista de violaciones, segun la severidad del criterio
    incumplido mas grave. Devuelve (etiqueta, [ids violados]).

    Jerarquia: numerico > critico > topologico > tecnologico/advertencia.
    Un punto puede ser a la vez DEGENERADO e INVIABLE; la etiqueta reporta la
    severidad dominante, pero la lista conserva los dos hallazgos, que es lo que
    exige la nota de redaccion del V3 (§10.6).
    """
    viol = [(i, s) for i, ok, s, _ in crit if not ok]
    sev = {s for _, s in viol}
    if 'numerico' in sev:
        etiqueta = 'NO_CONVERGE'
    elif 'critico' in sev:
        etiqueta = 'INVIABLE'
    elif 'topologico' in sev:
        etiqueta = 'DEGENERADO'
    elif 'tecnologico' in sev or 'advertencia' in sev:
        etiqueta = 'VALIDO con advertencia'
    else:
        etiqueta = 'VALIDO'
    return etiqueta, [i for i, _ in viol]


# ---------------------------------------------------------------------------
# SUPUESTOS ACTIVOS E INCERTIDUMBRE DECLARADA
#
# Regla del proyecto: ninguna variable de precision es obligatoria. Si falta,
# resolver() asume una hipotesis y esta funcion deja escrito que se asumio y que
# incertidumbre deja. Las cifras de efecto vienen de estudios ya medidos
# (resultados/incertidumbre-consolidada.md, presion_burbuja_o2_20260912.csv);
# no se recalculan en cada corrida.
# ---------------------------------------------------------------------------
def diagnostico_condensador(est, ind, par):
    """C_cold minimo (kW/K por kg/s de m_b) y caudal de agua minimo que exige
    la solucion ya resuelta. Una pasada de N evaluaciones de estado(). Devuelve
    dict(C_cold_min, m_agua_min) con None si el motor falla en el perfil."""
    Tk = par.get('T_agua_in', par['T_amb'])
    cp_w = par.get('cp_agua', 4.18)
    try:
        Ts, hs, n_salt = perfil_condensador(par['P_baja'], par['w_b'], est[9]['T'],
                                            est[9]['h'], est[8]['h'])
        cmin = C_cold_minimo(Ts, hs, Tk, par.get('dT_pp_cond', 0.0))
    except Exception:
        return dict(C_cold_min=None, m_agua_min=None, nodos_saltados=None)
    return dict(C_cold_min=cmin, m_agua_min=cmin / cp_w, nodos_saltados=n_salt)


def supuestos(est, ind, par, calcular_condensador=True):
    """Lista de dicts(variable, declarada, valor, efecto) con cada variable
    opcional: si se declaro, su valor; si no, la hipotesis asumida y la
    incertidumbre que deja. Incluye siempre la banda de composicion G4-01."""
    s = []
    add = lambda v, d, val, ef: s.append(dict(variable=v, declarada=d,
                                              valor=val, efecto=ef))
    deg = ind['degenerado']
    C_g = par.get('C_g')
    if C_g is not None and np.isfinite(C_g):
        add('C_g', True, f"{C_g:.4f} kW/K por kg/s de m_b", "fuente con capacidad finita")
    else:
        cmin = _cg_min_hrvg(est, ind, par)
        add('C_g', False, "fuente isoterma (H10)",
            "eta: sesgo < 0.3 % (Estudio 1; 0.0000 pp medido en la region valida "
            "mientras el pinzamiento no ata). Potencias totales, eta_recurso y "
            "T_gas,out no determinables. El regimen Kalina/degenerado puede "
            "cambiar con C_g finito."
            + (f" Minimo para no cruzar en el HRVG: C_g,min = {cmin:.3f} kW/K "
               f"por kg/s de m_b." if cmin is not None else ""))
    if par.get('dT_pp_gas') is None:
        add('dT_pp_gas', False, "0 K",
            "solo se exige no cruzar (limite termodinamico, no margen de diseno)")
    else:
        add('dT_pp_gas', True, f"{par['dT_pp_gas']:.2f} K", "margen de pinzamiento del HRVG")
    if par.get('T_min_gas') is None:
        add('T_min_gas', False, "no evaluado",
            "O6 (rocio acido) sin evaluar: no se puede fijar m_b,max ni la potencia total")
    else:
        add('T_min_gas', True, f"{par['T_min_gas']:.2f} K",
            "fija m_b,max; los MW son proporcionales a (T_f - T_min_gas)")

    C_cold = par.get('C_cold')
    if calcular_condensador:
        dc = diagnostico_condensador(est, ind, par)
    else:
        dc = dict(C_cold_min=None, m_agua_min=None)
    if dc['C_cold_min'] is None:
        txt_min = " (C_cold,min no calculado)"
    elif not np.isfinite(dc['C_cold_min']):
        txt_min = (" El extremo frio ya incumple el margen: ningun caudal de agua "
                   "alcanza el cierre declarado.")
    else:
        txt_min = (f" Exige C_cold >= {dc['C_cold_min']:.3f} kW/K por kg/s de m_b "
                   f"(~{dc['m_agua_min']:.2f} kg/s de agua por kg/s de m_b).")
    if C_cold is not None and np.isfinite(C_cold):
        add('C_cold', True, f"{C_cold:.4f} kW/K por kg/s de m_b",
            f"agua con capacidad finita; sale a {ind['T_agua_out']:.2f} K." + txt_min)
    else:
        add('C_cold', False, "sumidero infinito (el agua no se calienta)",
            "el estado 9 declarado solo es alcanzable con agua suficiente; con "
            "menos, T9 sube y O2 puede dejar de cumplirse." + txt_min)
    if 'T_agua_in' in par:
        add('T_agua_in', True, f"{par['T_agua_in']:.2f} K", "temperatura de entrada del agua")
    else:
        add('T_agua_in', False, f"T_amb = {par['T_amb']:.2f} K",
            "cada +1 K en la entrada del agua sube la P_baja minima que condensa "
            "~3-4 % (presion_burbuja_o2_20260912.csv, x_b 0.40-0.70)")
    if par.get('dT_pp_cond') is None:
        add('dT_pp_cond', False, "0 K",
            "solo se exige no cruzar en el condensador (limite termodinamico)")
    else:
        add('dT_pp_cond', True, f"{par['dT_pp_cond']:.2f} K", "margen interno del condensador")
    if 'T_agua_out_max' in par:
        add('T_agua_out_max', True, f"{par['T_agua_out_max']:.2f} K", "limite de vertido (O7)")
    else:
        add('T_agua_out_max', False, "no evaluado", "O7 (salida del agua) sin evaluar")

    if deg:
        comp = "nula en la rama degenerada (no hay composiciones de equilibrio)"
    else:
        comp = "2-4 % en eta en regimen Kalina (cota superior, Estudio 2)"
        if abs(ind['q4'] - 0.90) < 0.01:
            comp += (f". q4 = {ind['q4']:.4f} esta a menos de 0.01 del umbral de O1: "
                     "la etiqueta de O1 no es concluyente")
    add('composicion_G4-01', True, "banda +-0.01 de IAPWS G4-01 (siempre presente)", comp)
    return s


def imprimir_supuestos(sup):
    print("SUPUESTOS E INCERTIDUMBRE DECLARADA")
    for d in sup:
        marca = 'dado  ' if d['declarada'] else 'SUPUES'
        print(f"  [{marca}] {d['variable']:<17} {d['valor']}")
        print(f"           -> {d['efecto']}")


# ---------------------------------------------------------------------------
# DIMENSIONAMIENTO A LA FUENTE
#
# m_b no es un grado de libertad termodinamico (todas las ecuaciones del ciclo
# son especificas; m_b solo multiplica Q_i, W_t, W_p, W_net). Lo unico que
# importa es la razon C_g/m_b. "Dimensionar m_b a la fuente" es elegir esa
# razon: m_b sale reportado como el tamano de planta resultante.
# ---------------------------------------------------------------------------
def dimensionar(par, C_g_absoluto, T_min_gas, dT_pp_gas=0.0,
                 mb_lo=0.05, mb_hi=None, C_cold_absoluto=None):
    """Mayor m_b (caudal de solucion basica, kg/s) que cumple a la vez:

      (1) T_gas,out >= T_min_gas          (rocio acido + margen de diseno)
      (2) DT >= dT_pp_gas en todo el HRVG SIN que el limitador de pinzamiento
          tenga que recortar h2 -- el mayor m_b para el que el cierre
          DECLARADO (dT_app o eps_HRVG) sigue describiendo el equipo.

    Manda el que ate primero (el m_b mas chico de los dos). C_g_absoluto es la
    capacidad calorifica REAL de la fuente [kW/K], un numero fijo del equipo;
    para evaluar el ciclo con un m_b de prueba se normaliza como
    C_g_absoluto/m_b, porque resolver() trabaja en caudales especificos
    (normalizados a m_b=1 -- ver docstring del modulo).

    Devuelve dict(m_b, cual_ato, m_b_rocio, m_b_pinzamiento, est, ind).

    Nota: la busqueda de raiz es una bisatanteo simple con expansion de
    intervalo; para C_g_absoluto o T_min_gas en regimenes extremos (sin
    solucion factible en ningun m_b, o factible en todo el rango explorado)
    devuelve el extremo del intervalo en vez de fallar -- revisar `cual_ato`
    y los valores de m_b_rocio/m_b_pinzamiento antes de usar el resultado.
    """
    base = dict(par)

    def evalua(mb):
        p = dict(base)
        p['C_g'] = C_g_absoluto / mb
        p['dT_pp_gas'] = dT_pp_gas
        if C_cold_absoluto is not None:     # agua real de la planta [kW/K]
            p['C_cold'] = C_cold_absoluto / mb
        return resolver(p)

    def g_rocio(mb):
        _, ind = evalua(mb)
        return ind['T_gas_out'] - T_min_gas

    def g_pinch(mb):
        _, ind = evalua(mb)
        return ind['DT_min_HRVG'] - dT_pp_gas

    def raiz(g, lo, hi_ini):
        hi = hi_ini
        flo, fhi = g(lo), g(hi)
        n = 0
        while flo * fhi > 0 and n < 30:
            hi *= 1.5
            fhi = g(hi)
            n += 1
        if flo * fhi > 0:
            return hi if flo > 0 else lo   # no ata en el rango explorado
        return brentq(g, lo, hi, xtol=1e-4)

    if mb_hi is None:
        mb_hi = max(C_g_absoluto, 1.0)
    mb_r = raiz(g_rocio, mb_lo, mb_hi)
    mb_p = raiz(g_pinch, mb_lo, mb_hi)
    mb = min(mb_r, mb_p)
    cual = 'rocio_acido' if mb_r <= mb_p else 'pinzamiento'
    est, ind = evalua(mb)
    return dict(m_b=mb, cual_ato=cual, m_b_rocio=mb_r, m_b_pinzamiento=mb_p,
                est=est, ind=ind)
