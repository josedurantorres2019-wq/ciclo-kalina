"""PRUEBA 4 — reproduccion de un punto publicado de Elsayed et al. (2013)
por el METODO A (approach + pinch).

Punto: eta = 11.38 % a 15 bar, x_NH3 = 0.55 masica, fuente 373 K, sumidero 283 K
       (texto de la seccion 3 del paper, sobre la Figura 4).
Anclas independientes publicadas en el mismo parrafo: eta_Carnot = 24 %,
eta_II = 47 %.

Hipotesis suplidas (el paper no las da) — congeladas ANTES de ver resultados:
  H-a  approach del HRVG = 4 K   -> T2 = 369 K
  H-b  estado 9 = liquido saturado a T_amb + 4 K = 287 K -> P_baja = P_burbuja
  H-c  pinch del regenerador de 4 K cerrado en el extremo frio (T6 = T10 + 4 K),
       con verificacion posterior del extremo caliente
  H-d  bomba por via entropica en lugar de v*dP (desviacion medida: 0.05 %)

Criterio de aceptacion declarado de antemano: |eta_calc - 11.38| <= 1.5 puntos.
"""
import numpy as np
import nh3h2o as ng
import kalina as k

W_B = 0.55
T_F = 373.0
T_AMB = 283.0
P_ALTA = 1.5            # MPa
DT_APP = 4.0
DT_PP = 4.0
DT_COND = 4.0
ETA_T = ETA_P = 0.80
ETA_PUB = 11.38         # %
TOL_PUB = 1.5           # puntos porcentuales


def caso(dT_app=DT_APP, dT_cond=DT_COND, P_baja=None, verbose=True):
    xb = ng.w2m(W_B)
    T9 = T_AMB + dT_cond
    if P_baja is None:
        P_baja = ng.bubbleP(T9, xb)[0]
    par = dict(P_alta=P_ALTA, P_baja=P_baja, w_b=W_B, T_f=T_F, T_amb=T_AMB,
               eta_t=ETA_T, eta_p=ETA_P, metodo='A',
               cierre_HRVG=dict(dT_app=dT_app),
               cierre_reg=dict(dT_pp=DT_PP),
               cierre_cond=dict(dT_app=dT_cond))
    est, ind = k.resolver(par)
    crit = k.criterios(est, ind, par)
    if verbose:
        reporte(par, est, ind, crit)
    return par, est, ind, crit


def reporte(par, est, ind, crit):
    print("=" * 78)
    print("PRUEBA 4 — Elsayed et al. (2013), 15 bar / x_b=0.55 / 373 K / 283 K")
    print("=" * 78)
    print(f"P_alta = {par['P_alta']*1000:.1f} kPa    "
          f"P_baja = {par['P_baja']*1000:.2f} kPa (resultado de H-b)")
    print(f"T_f = {par['T_f']:.2f} K   T_amb = {par['T_amb']:.2f} K   "
          f"eta_t = eta_p = {par['eta_t']:.2f}")
    print(f"x_b = {W_B:.4f} masica = {ng.w2m(W_B):.4f} molar")
    print()
    print("ESTADOS")
    print(f"{'#':>3} {'T [C]':>9} {'P [kPa]':>9} {'w_NH3':>7} {'h [kJ/kg]':>11} "
          f"{'s [kJ/kgK]':>11} {'q':>8}  fase")
    for i in range(1, 11):
        e = est[i]
        print(f"{i:>3} {e['T']-273.15:>9.3f} {e['P']*1000:>9.2f} {e['w']:>7.4f} "
              f"{e['h']:>11.3f} {e['s']:>11.5f} {e['q']:>8.4f}  {e['fase']}")
    print()
    print("CAUDALES Y COMPOSICIONES DEL SEPARADOR")
    print(f"  m_r/m_b = {ind['m_r']:.5f}   m_p/m_b = omega = {ind['m_p']:.5f}")
    print(f"  w_r = {ind['w_r']:.5f}   w_p = {ind['w_p']:.5f}   q2 = {ind['q2']:.5f}")
    print()
    print("PAR DE CIERRE — las cuatro columnas, siempre")
    print(f"  eps_HRVG = {ind['eps_HRVG']:.5f} (despejada)   "
          f"dT_app = {ind['dT_app']:.3f} K (declarada)")
    print(f"  eps_reg  = {ind['eps_reg']:.5f} (despejada)   "
          f"dT_pp,reg = {ind['dT_pp']:.3f} K (declarada)")
    print()
    print("INDICADORES  [kJ/kg de fluido basico; con m_b = 1 kg/s son kW]")
    print(f"  Q_i     = {ind['Qi']:10.3f}")
    print(f"  Q_out   = {ind['Qout']:10.3f}")
    print(f"  W_t     = {ind['Wt']:10.3f}")
    print(f"  W_p     = {ind['Wp']:10.3f}")
    print(f"  W_net   = {ind['Wnet']:10.3f}")
    print(f"  eta     = {100*ind['eta']:10.3f} %")
    print(f"  eta_Car = {100*ind['eta_carnot']:10.3f} %")
    print(f"  eta_II  = {100*ind['eta_II']:10.3f} %")
    print(f"  RBW     = {100*ind['RBW']:10.3f} %")
    print()
    print("CRITERIOS")
    for i, ok, sev, txt in crit:
        print(f"  [{'OK ' if ok else 'NO '}] {i:<4} {sev:<13} {txt}")
    print()
    print("CONTRASTE CON EL PAPER")
    d = 100 * ind['eta'] - ETA_PUB
    print(f"  eta publicada    = {ETA_PUB:.2f} %")
    print(f"  eta calculada    = {100*ind['eta']:.2f} %")
    print(f"  diferencia       = {d:+.2f} puntos   "
          f"(criterio declarado: |d| <= {TOL_PUB:.1f})")
    print(f"  VEREDICTO PRUEBA 4: {'REPRODUCE' if abs(d) <= TOL_PUB else 'NO REPRODUCE'}")
    print(f"  ancla independiente eta_Carnot: publicada 24 %, "
          f"calculada {100*ind['eta_carnot']:.2f} %")
    print(f"  ancla independiente eta_II    : publicada 47 %, "
          f"calculada {100*ind['eta_II']:.2f} %")


def sensibilidad_Pbaja(P_kPa=(200, 250, 273.5, 300, 350, 400, 500, 600)):
    print()
    print("=" * 78)
    print("SENSIBILIDAD A LA HIPOTESIS H-b (P_baja no declarada por el paper)")
    print("=" * 78)
    print(f"{'P_baja [kPa]':>13} {'T9 [C]':>9} {'eta [%]':>9} {'W_net':>9} "
          f"{'q4':>7} {'fase 9':>10}")
    for Pk in P_kPa:
        try:
            par, est, ind, _ = caso(P_baja=Pk / 1000.0, verbose=False)
            print(f"{Pk:>13.1f} {est[9]['T']-273.15:>9.3f} {100*ind['eta']:>9.3f} "
                  f"{ind['Wnet']:>9.3f} {ind['q4']:>7.4f} {est[9]['fase']:>10}")
        except Exception as ex:
            print(f"{Pk:>13.1f}   fallo: {ex}")


def sensibilidad_approach(dTs=(0.0, 2.0, 4.0, 6.0, 8.0)):
    print()
    print("=" * 78)
    print("SENSIBILIDAD A LA HIPOTESIS H-a (approach del HRVG no declarado)")
    print("=" * 78)
    print(f"{'dT_app [K]':>11} {'T2 [C]':>9} {'q2':>8} {'omega':>8} "
          f"{'eta [%]':>9} {'W_net':>9}")
    for d in dTs:
        try:
            par, est, ind, _ = caso(dT_app=d, verbose=False)
            print(f"{d:>11.1f} {est[2]['T']-273.15:>9.3f} {ind['q2']:>8.4f} "
                  f"{ind['omega']:>8.4f} {100*ind['eta']:>9.3f} {ind['Wnet']:>9.3f}")
        except Exception as ex:
            print(f"{d:>11.1f}   fallo: {ex}")


if __name__ == "__main__":
    caso()
    sensibilidad_Pbaja()
    sensibilidad_approach()
