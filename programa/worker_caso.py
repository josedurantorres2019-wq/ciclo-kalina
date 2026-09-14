"""Script worker que resuelve UN caso del ciclo Kalina en un proceso aparte.

Se invoca como subprocess (no via multiprocessing.Process) desde motor_ui.py:
en Windows, multiprocessing con metodo 'spawn' reimporta el script principal
para reconstruir el entorno del proceso hijo, y ese script principal aqui es
la app de Streamlit -- reimportarla en el hijo es fragil e innecesariamente
costoso. Un subprocess.run() normal con `python worker_caso.py` no tiene ese
problema: arranca limpio, sin nada de Streamlit cargado.

Contrato: lee un JSON por stdin con {"par": {...}, "calcular_condensador": bool},
escribe un JSON por stdout con el resultado. No imprime nada mas por stdout.
"""
import json
import sys
import traceback
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def _limpiar_valor(v):
    """Convierte tipos numpy a tipos nativos de Python serializables en JSON."""
    tipo = type(v).__name__
    if tipo.startswith("float") and tipo != "float":
        return float(v)
    if tipo.startswith(("int", "uint")) and tipo not in ("int",):
        return int(v)
    if tipo == "bool_":
        return bool(v)
    return v


def _limpiar_dict(d):
    return {k: _limpiar_valor(v) for k, v in d.items()}


def main():
    entrada = json.loads(sys.stdin.read())
    par = entrada["par"]
    calcular_condensador = entrada.get("calcular_condensador", False)
    try:
        import kalina as k
        est, ind = k.resolver(par)
        crit = k.criterios(est, ind, par)
        etiqueta, violados = k.clasificar(crit)
        sup = k.supuestos(est, ind, par, calcular_condensador=calcular_condensador)
        salida = dict(
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
        salida = dict(estado=tipo, mensaje=str(e), traceback=traceback.format_exc())
    sys.stdout.write(json.dumps(salida))


if __name__ == "__main__":
    main()
