"""ACOTACION DE LA INCERTIDUMBRE DE LOS VALORES ABSOLUTOS.

  ESTUDIO 1 — capacidad calorifica finita de la fuente.
    Convierte el sesgo no cuantificado de la hipotesis H10 (fuente como
    reservorio) en un intervalo declarado y medido, usando el limitador de
    pinzamiento de kalina.py (verificacion de no cruce interno de
    temperaturas en el HRVG, que el V3 dejaba como extension).

  ESTUDIO 2 — banda de incertidumbre de la composicion de equilibrio.
    Propaga los +-0.01 que la guia IAPWS G4-01 declara para la composicion de
    las fases en equilibrio y mide su efecto sobre m_r, W_net y eta.

--- ACTUALIZACION 2026-09-10: reescrito sobre kalina.resolver(par['C_g']) ---

Este archivo tenia su PROPIA implementacion del limitador de pinzamiento
(perfil/C_g_minimo/h2_maximo, mas un lazo externo de hasta 3 pasadas en
estudio_1 para converger T1). Se elimino: kalina.py ahora hace exactamente
eso -y mejor- DENTRO de resolver(), resolviendo h1 y h2 juntos en el mismo
lazo en vez de por fuera con iteraciones manuales. Declarar par['C_g'] (y
opcionalmente par['dT_pp_gas']) es todo lo que hace falta; el resto de este
docstring describe la fisica, que no cambio.

**Correccion C26 -- direccion del gas.** La formula de esta version anterior
(`T_g(h) = T_f - (h-h1)/C_g`, gas entrando caliente por el extremo T1) tenia
el flujo en PARALELO (co-corriente), inconsistente con el cierre por approach
`dT_app = T_f - T2` que ya supone contracorriente (T2, no T1, es el que se
acerca a T_f). Corregido en kalina.pinzamiento(): `T_g(h) = T_f - (h2-h)/C_g`.
Esto INVALIDA los numeros que este estudio reportaba antes de la correccion
(el regimen volteando a Kalina a ~10-15 kg/s de gas) -- ver
../../resultados/hallazgo-degeneracion-caso-base.md para el detalle y lo que
queda por verificar con la formula corregida.

--- POR QUE EL PLANTEAMIENTO eps-NTU NO BASTA (correccion al primer intento) ---

La primera version acotaba el efecto de C_g por la via eps-NTU: Q_max =
C_min*(T_h,ent - T_c,ent), de modo que al bajar C_g por debajo de C_wf el calor
maximo transferible cae y con el la eficiencia. Es correcto pero INOPERANTE
aqui: para phi = C_g/C_wf >= 1 el C_min sigue siendo el fluido de trabajo, Q_max
no cambia, y eta sale IDENTICA para phi = 1, 2, 10 o infinito. El sesgo que
queremos acotar no aparece por ningun lado.

La razon fisica es que en una mezcla zeotropica la restriccion que muerde no es
la capacidad TOTAL del gas sino el PERFIL: la curva h(T) del fluido de trabajo
es fuertemente no lineal (tramo de liquido, evaporacion con deslizamiento de
86 K, tramo de vapor) mientras que la del gas es una recta. Dos curvas de
distinta forma pueden tener area suficiente y aun asi tocarse o cruzarse en el
interior. Ese cruce es una violacion de la segunda ley que S5 (T2 <= T_f) no ve,
exactamente como advierte la correccion C12 del V3.

Planteamiento adoptado, que si mide el sesgo (implementado en
kalina.perfil_HRVG / kalina.pinzamiento, dentro de resolver()):

  1. Se discretiza el HRVG en N nodos EN TEMPERATURA entre T1 y T2 (barato: no
     hay que invertir h->T, se evalua h(T) directamente) y se obtiene el perfil
     real del fluido de trabajo.
  2. El gas es lineal en el calor transferido y entra caliente por el extremo
     T2 (contracorriente, correccion C26):  T_g(h) = T_f - (h2 - h)/C_g.
  3. Para cada C_g se busca el h2 MAXIMO compatible con  DT(h) >= DT_pp  en todo
     el equipo, y se recorta el cierre (approach o efectividad, el que se haya
     declarado) a ese h2 si hace falta. Ya no importa que metodo se declare:
     el limitador actua DESPUES del cierre, sobre el resultado de cualquiera.
  4. h1 y h2 se resuelven JUNTOS en el mismo lazo exterior de resolver() (antes
     era un lazo externo de hasta 3 pasadas, en este archivo).
  5. La diferencia contra el caso C_g=None (fuente como reservorio) ES el
     sesgo de H10, medido.

Como subproducto se obtiene C_g,min (`ind['Cg_min_HRVG']`, informativo, solo
cuando C_g no se declara): la capacidad calorifica minima de la fuente que
hace el intercambiador termodinamicamente posible con el Q_i que el modelo
actual predice. Por debajo de ese valor, el resultado con fuente isoterma no
solo es impreciso: es imposible.

--- ANCLAJE DEL CONTEXTO ---

Fuente: gases de escape industriales a 350 C. Se descarta el anclaje
geotermico: los recursos geotermicos convencionales estan entre 100 y 250 C y
el portador es salmuera liquida (c_p ~ 4.2 kJ/kg-K), no gas; solo pozos
supercriticos experimentales alcanzan 350 C.

c_p de gases de combustion a 300-400 C: 1.10 kJ/kg-K (tablas de productos de
combustion / aire, Cengel y Boles, 2012). Se usa SOLO para traducir C_g a un
caudal de gas y juzgar plausibilidad; el analisis esta hecho en C_g, que es la
variable que gobierna el problema.

DT_pp del lado gas: 10 K como valor de diseno y 0 K como limite termodinamico.
El rango habitual de aproximacion minima gas-liquido en recuperacion de calor
es de 8 a 15 K (practica de diseno de calderas de recuperacion); por debajo de
eso el area requerida crece sin cota.

La normalizacion es por kg/s de fluido de trabajo (m_b = 1 kg/s), de modo que
C_g y m_gas se leen como "por cada kg/s de mezcla amoniaco-agua".
"""
import numpy as np
import nh3h2o as ng
import kalina as k

CP_GAS = 1.10           # kJ/kg-K, gases de combustion a 300-400 C
BANDA_W = 0.01          # banda de composicion declarada por IAPWS G4-01
DT_PP_GAS = 10.0        # K, aproximacion minima de diseno en el lado gas

BASE = dict(P_alta=3.0, P_baja=0.4, w_b=0.50, T_f=623.15, T_amb=300.033,
            eta_t=0.85, eta_p=0.75, metodo='B',
            cierre_HRVG=dict(eps=0.85), cierre_reg=dict(eps=0.75),
            cierre_cond=dict(eps=0.80))


def caso(base=None, perturb=(0.0, 0.0), T2=None):
    par = {kk: (dict(vv) if isinstance(vv, dict) else vv)
           for kk, vv in (base or BASE).items()}
    if T2 is not None:                      # cierre del HRVG por approach
        par['cierre_HRVG'] = dict(dT_app=par['T_f'] - T2)
    par['perturb_w'] = perturb
    est, ind = k.resolver(par)
    return par, est, ind


# ---------------------------------------------------------------------------
# ESTUDIO 1
# ---------------------------------------------------------------------------
def estudio_1(m_gas=(80.0, 40.0, 25.0, 15.0, 10.0), dT_pp=DT_PP_GAS, base=None):
    """C_g se declara directamente a resolver() (par['C_g'], par['dT_pp_gas']):
    h1 y h2 se resuelven juntos dentro del lazo exterior de kalina.resolver(),
    con el limitador de pinzamiento (contracorriente, correccion C26) actuando
    despues del cierre que se haya declarado en `base` (approach o
    efectividad). Barrido con ARRANQUE TIBIO: cada punto siembra el h1 del
    punto anterior (par['h1_semilla']) -mismo principio que estado_de()-, que
    en la practica da ~15-20x de mejora salvo justo en el entorno de un cruce
    de fase, donde la convergencia sigue siendo dificil (ver
    ../../resultados/hallazgo-degeneracion-caso-base.md).
    """
    print("=" * 100)
    print("ESTUDIO 1 — capacidad calorifica finita de la fuente (hipotesis H10)")
    print("=" * 100)
    par0, est0, ind0 = caso(base=base)
    Tf = par0['T_f']

    print(f"\nREFERENCIA (C_g no declarado, hipotesis H10):")
    print(f"  regimen {'DEGENERADO' if ind0['degenerado'] else 'KALINA'}   "
          f"T1 = {est0[1]['T']-273.15:.2f} C   T2 = {est0[2]['T']-273.15:.2f} C   "
          f"q2 = {ind0['q2']:.4f}")
    print(f"  Q_i = {ind0['Qi']:.2f}   W_net = {ind0['Wnet']:.2f}   "
          f"eta = {100*ind0['eta']:.3f} %")
    Cg_min_dis = ind0['Cg_min_HRVG']
    print(f"\nCAPACIDAD MINIMA DE LA FUENTE para sostener ese Q_i sin cruce interno:")
    print(f"  C_g,min = {Cg_min_dis:.2f} kW/K  ->  m_gas = {Cg_min_dis/CP_GAS:.1f} "
          f"kg/s      (DT_pp = {dT_pp:.0f} K, criterio de diseno)")
    print(f"\nCICLO CORREGIDO para cada fuente real (T2 limitado por el "
          f"pinzamiento interno, DT_pp = {dT_pp:.0f} K):")
    print(f"{'m_gas':>7} {'C_g':>8} {'T2 [C]':>9} {'q2':>8} {'regimen':>11} "
          f"{'Q_i':>9} {'W_net':>8} {'eta %':>8} {'d eta %':>9} "
          f"{'Tg,out C':>9} {'@z':>5} {'recorte':>7} {'et_rec %':>9}")

    filas = []
    h1_seed = ind0['h1']
    for mg in m_gas:
        C_g = mg * CP_GAS
        par = {kk: (dict(vv) if isinstance(vv, dict) else vv)
               for kk, vv in (base or BASE).items()}
        par['C_g'] = C_g
        par['dT_pp_gas'] = dT_pp
        par['h1_semilla'] = h1_seed
        try:
            est, ind = k.resolver(par)
            h1_seed = ind['h1']
            d = 100 * (ind['eta'] - ind0['eta']) / ind0['eta']
            z = ind['z_min_HRVG']
            print(f"{mg:>7.1f} {C_g:>8.2f} {est[2]['T']-273.15:>9.2f} "
                  f"{ind['q2']:>8.4f} "
                  f"{('DEGENERADO' if ind['degenerado'] else 'KALINA'):>11} "
                  f"{ind['Qi']:>9.2f} {ind['Wnet']:>8.2f} {100*ind['eta']:>8.3f} "
                  f"{d:>+9.2f} {ind['T_gas_out']-273.15:>9.2f} "
                  f"{(z if z is not None else float('nan')):>5.2f} "
                  f"{str(ind['pinch_recortado']):>7} {100*ind['eta_recurso']:>9.3f}")
            filas.append((mg, C_g, est, ind))
        except Exception as ex:
            print(f"{mg:>7.1f} {C_g:>8.2f}   {type(ex).__name__}: {ex}")
    return ind0, filas


def tabla_perfil(base=None, N=15):
    """Perfil interno del HRVG (caso de referencia, C_g no declarado), para ver
    DONDE esta el codo de la curva. Usa kalina.perfil_HRVG directamente."""
    par0, est0, ind0 = caso(base=base)
    Tf, Ph, wb = par0['T_f'], par0['P_alta'], par0['w_b']
    Ts, hs = k.perfil_HRVG(Ph, wb, est0[1]['T'], est0[1]['h'],
                            est0[2]['T'], est0[2]['h'], N)
    print()
    print("PERFIL INTERNO DEL FLUIDO DE TRABAJO EN EL HRVG (caso de referencia)")
    print("(contracorriente: el gas entraria caliente por el extremo T2 -- correccion C26)")
    print(f"{'z=Q/Qi':>8} {'T_wf [C]':>10} {'h [kJ/kg]':>11} "
          f"{'C_g req (DT=0)':>15} {'fase':>10}")
    h2 = hs[-1]
    for T, h in zip(Ts, hs):
        z = (h - hs[0]) / (hs[-1] - hs[0])
        req = (h2 - h) / (Tf - T) if T < Tf else np.inf
        e = k.estado(float(T), Ph, wb)
        print(f"{z:>8.3f} {T-273.15:>10.2f} {h:>11.2f} {req:>15.2f} {e['fase']:>10}")


# ---------------------------------------------------------------------------
# ESTUDIO 2
# ---------------------------------------------------------------------------
def estudio_2(base=None, T2=None):
    print()
    print("=" * 100)
    print("ESTUDIO 2 — banda de composicion de equilibrio de IAPWS G4-01 (+-0.01)")
    print("=" * 100)
    par0, est0, ind0 = caso(base=base, T2=T2)
    if ind0['degenerado']:
        print("  Caso DEGENERADO: no hay separacion, w_r = w_p = w_b, y la banda de")
        print("  composicion no tiene sobre que actuar. Su contribucion a la")
        print("  incertidumbre es EXACTAMENTE CERO en este regimen.")
        return 0.0
    print(f"{'dw_r':>7} {'dw_p':>7} {'w_r':>8} {'w_p':>8} {'m_r/m_b':>9} "
          f"{'W_net':>9} {'eta %':>8} {'d eta %':>9} {'res.sep':>10}")
    combos = [(0, 0), (+BANDA_W, 0), (-BANDA_W, 0), (0, +BANDA_W), (0, -BANDA_W),
              (+BANDA_W, -BANDA_W), (-BANDA_W, +BANDA_W)]
    ref, peor, peor_mr = None, 0.0, 0.0
    for dwr, dwp in combos:
        par, est, ind = caso(base=base, perturb=(dwr, dwp), T2=T2)
        if ref is None:
            ref, ref_mr = ind['eta'], ind['m_r']
        d = 100 * (ind['eta'] - ref) / ref
        peor = max(peor, abs(d))
        peor_mr = max(peor_mr, abs(100 * (ind['m_r'] - ref_mr) / ref_mr))
        print(f"{dwr:>+7.3f} {dwp:>+7.3f} {ind['w_r']:>8.4f} {ind['w_p']:>8.4f} "
              f"{ind['m_r']:>9.5f} {ind['Wnet']:>9.3f} {100*ind['eta']:>8.4f} "
              f"{d:>+9.3f} {ind['res_separador']:>10.3e}")
    print()
    print(f"  Desviacion maxima dentro de la banda:  eta {peor:.3f} %   "
          f"m_r/m_b {peor_mr:.3f} %")
    return peor


if __name__ == "__main__":
    estudio_1()
    tabla_perfil()
    estudio_2()
