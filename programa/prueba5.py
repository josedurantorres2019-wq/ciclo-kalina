"""PRUEBA 5 — consistencia cruzada A <-> B anclada en el punto de Elsayed.

Que hace:
  1. Resuelve el punto de la prueba 4 por el METODO A (approach + pinch).
  2. Lee del resultado las tres efectividades EQUIVALENTES (HRVG, regenerador,
     condensador) — la conversion biyectiva de la regla del grado de libertad.
  3. Resuelve el MISMO punto por el METODO B con esas efectividades, que obliga
     al resolvedor a recorrer el lazo exterior sobre h1 (Brent) y el lazo
     interior frio, ninguno de los cuales se ejerce en el metodo A.
  4. Compara estado por estado.

Que certifica: que los dos cierres son la misma fisica y que el resolvedor de
lazo es correcto. Como el punto es el mismo que valido la prueba 4 contra el
paper, la validacion externa del metodo A se TRANSFIERE al metodo B.
Que NO certifica: nada sobre el rango de operacion del proyecto (350 C).

Control de poder discriminante (paso 5): se repite la comparacion con una
efectividad perturbada. Si la prueba no fallara ahi, no estaria midiendo nada.
"""
import numpy as np
import nh3h2o as ng
import kalina as k
import prueba4 as p4

TOL_T = 1e-3        # K
TOL_H = 1e-3        # kJ/kg
TOL_REL = 1e-6      # indicadores, relativo


def par_B(par_A, ind_A, d_eps_HRVG=0.0, d_eps_reg=0.0, d_eps_cond=0.0):
    """Traduce la corrida A a su equivalente en el metodo B."""
    q = dict(par_A)
    q['metodo'] = 'B'
    q['cierre_HRVG'] = dict(eps=ind_A['eps_HRVG'] + d_eps_HRVG)
    q['cierre_reg'] = dict(eps=ind_A['eps_reg'] + d_eps_reg)
    q['cierre_cond'] = dict(eps=ind_A['eps_cond'] + d_eps_cond)
    return q


def comparar(estA, indA, estB, indB, titulo, tol_T=TOL_T, tol_h=TOL_H,
             tol_rel=TOL_REL):
    print("=" * 78)
    print(titulo)
    print("=" * 78)
    print(f"{'#':>3} {'T_A [K]':>11} {'T_B [K]':>11} {'dT [K]':>11} "
          f"{'h_A':>11} {'h_B':>11} {'dh':>11}")
    ok = True
    for i in range(1, 11):
        a, b = estA[i], estB[i]
        dT, dh = b['T'] - a['T'], b['h'] - a['h']
        ok &= abs(dT) <= tol_T and abs(dh) <= tol_h
        print(f"{i:>3} {a['T']:>11.5f} {b['T']:>11.5f} {dT:>11.2e} "
              f"{a['h']:>11.4f} {b['h']:>11.4f} {dh:>11.2e}")
    print()
    print(f"{'indicador':>12} {'metodo A':>14} {'metodo B':>14} {'dif. rel.':>12}")
    for key, nom in (('Qi', 'Q_i'), ('Qout', 'Q_out'), ('Wt', 'W_t'),
                     ('Wp', 'W_p'), ('Wnet', 'W_net'), ('eta', 'eta'),
                     ('m_r', 'm_r/m_b'), ('w_r', 'w_r'), ('w_p', 'w_p'),
                     ('q2', 'q2'), ('q4', 'q4')):
        a, b = indA[key], indB[key]
        rel = abs(b - a) / max(abs(a), 1e-12)
        ok &= rel <= tol_rel
        print(f"{nom:>12} {a:>14.7f} {b:>14.7f} {rel:>12.2e}")
    print()
    print("Par de cierre — cada metodo declara uno y despeja el otro:")
    print(f"{'':>14} {'declarado':>14} {'despejado (A)':>16} {'despejado (B)':>16}")
    print(f"{'HRVG':>14} {'dT_app = 4 K':>14} "
          f"{indA['eps_HRVG']:>16.7f} {indB['dT_app']:>13.5f} K")
    print(f"{'regenerador':>14} {'dT_pp = 4 K':>14} "
          f"{indA['eps_reg']:>16.7f} {indB['dT_pp']:>13.5f} K")
    print(f"{'condensador':>14} {'dT_app = 4 K':>14} "
          f"{indA['eps_cond']:>16.7f} {indB['dT_app_cond']:>13.5f} K")
    print()
    print(f"lazo exterior (metodo B): ancho final del intervalo de Brent = "
          f"{indB['ancho_h1']:.1e} kJ/kg   [criterio N1]")
    print(f"lazo interior frio      : {indB['iter_frio']} pasadas "
          f"(metodo A: {indA['iter_frio']})")
    print()
    print(f"VEREDICTO: {'COINCIDEN' if ok else 'NO COINCIDEN'}")
    return ok


def main():
    print("### Paso 1-2: metodo A y lectura de las efectividades equivalentes\n")
    parA, estA, indA, critA = p4.caso(verbose=False)
    print(f"  eps_HRVG = {indA['eps_HRVG']:.10f}")
    print(f"  eps_reg  = {indA['eps_reg']:.10f}")
    print(f"  eps_cond = {indA['eps_cond']:.10f}")
    print()

    print("### Paso 3-4: metodo B con esas efectividades\n")
    parB = par_B(parA, indA)
    estB, indB = k.resolver(parB)
    ok = comparar(estA, indA, estB, indB,
                  "PRUEBA 5 — consistencia cruzada A <-> B (anclada)")

    print()
    print("### Criterios evaluados sobre la solucion del metodo B\n")
    for i, o, sev, txt in k.criterios(estB, indB, parB):
        print(f"  [{'OK ' if o else 'NO '}] {i:<4} {sev:<13} {txt}")

    print()
    print("### Paso 5: control de poder discriminante\n")
    print("Se repite la comparacion con eps_HRVG perturbada en +0.001.")
    print("Si el test siguiera dando COINCIDEN, no estaria midiendo nada.\n")
    parB2 = par_B(parA, indA, d_eps_HRVG=1e-3)
    estB2, indB2 = k.resolver(parB2)
    comparar(estA, indA, estB2, indB2,
             "CONTROL — metodo B con eps_HRVG + 0.001 (debe FALLAR)")
    return ok


if __name__ == "__main__":
    main()
