"""Script worker que resuelve una LISTA de casos del ciclo Kalina en UN SOLO
proceso persistente (a diferencia de worker_caso.py, que resuelve uno solo).

Por que hace falta: kalina.py mantiene cache a nivel de modulo (_puros,
_flash, _env, _ultimaT) que se calienta con cada llamada y acelera las
siguientes -- es la razon por la que los barridos del propio proyecto
(barrido_sensibilidad.py, etc.) corren muchos puntos en un solo proceso y
solo el primero es lento. Resolver cada caso de la rejilla en su PROPIO
subprocess (como hace worker_caso.py) tira ese cache a la basura en cada
caso: se midio ~20-30 s por caso incluso en la rama degenerada barata, contra
un puñado de segundos para el primero y fracciones de segundo para el resto
dentro de un mismo proceso. motor_ui.resolver_grid() usa este script para
correr toda la rejilla (o el tramo que quede) en un proceso compartido, y
solo lo reinicia (perdiendo el cache) si un caso concreto se cuelga mas alla
del timeout -- ver esa funcion para la logica de reinicio.

Contrato: lee por stdin un JSON con una lista de
[{"par": {...}, "calcular_condensador": bool}, ...]. Por cada caso, EN ORDEN,
escribe una linea JSON por stdout con el resultado (mismo formato por caso
que worker_caso.py) y hace flush inmediato, para que el proceso padre pueda
leer el progreso incrementalmente y aplicar un timeout por caso sin esperar a
que termine el lote completo.
"""
import json
import sys
import traceback
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from worker_caso import _limpiar_dict  # mismo directorio


def _resolver_uno(par, calcular_condensador):
    try:
        import kalina as k
        est, ind = k.resolver(par)
        crit = k.criterios(est, ind, par)
        etiqueta, violados = k.clasificar(crit)
        sup = k.supuestos(est, ind, par, calcular_condensador=calcular_condensador)
        return dict(
            estado="ok",
            est={str(i): _limpiar_dict(e) for i, e in est.items()},
            ind=_limpiar_dict(ind),
            crit=[[cid, bool(ok), sev, txt] for cid, ok, sev, txt in crit],
            etiqueta=etiqueta,
            violados=list(violados),
            sup=sup,
        )
    except Exception as e:
        tipo = "no_converge" if type(e).__name__ == "NoConverge" else "error"
        return dict(estado=tipo, mensaje=str(e), traceback=traceback.format_exc())


def main():
    casos = json.loads(sys.stdin.read())
    # Arranque tibio encadenado: kalina.resolver() acepta par['h1_semilla']
    # (opcional) para acotar el intervalo de Brent del lazo exterior alrededor
    # de un h1 vecino ya resuelto, en vez de partir del rango fisico completo
    # -- ver su docstring. Sin esto, CADA caso hace una busqueda de Brent a
    # ciegas sobre todo el rango, que medido puede tardar varios minutos por
    # punto incluso sin fuente de capacidad finita declarada (ver
    # plan-interfaz-web-kalina.md, actualizacion sobre resolver_grid). No
    # sustituye un h1_semilla que el propio caso ya traiga declarado.
    h1_semilla = None
    for caso in casos:
        par = dict(caso["par"])
        if h1_semilla is not None and "h1_semilla" not in par:
            par["h1_semilla"] = h1_semilla
        salida = _resolver_uno(par, caso.get("calcular_condensador", False))
        if salida["estado"] == "ok":
            h1_semilla = salida["ind"]["h1"]
        sys.stdout.write(json.dumps(salida) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
