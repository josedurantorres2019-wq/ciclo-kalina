"""Genera los fixtures de regresion numerica (golden master).

Corre los casos REALES actuales y congelo sus salidas en tests/fixtures/*.json.
Esos JSON son la verdad de referencia: antes de tocar el motor hay que
regenerarlos y revisar con git diff que NO cambian.

Uso (desde programa/):
    .venv\\Scripts\\python -m tests.generar_fixtures            # todo
    .venv\\Scripts\\python -m tests.generar_fixtures --solo motor  # solo estados+flashes
    .venv\\Scripts\\python -m tests.generar_fixtures --solo A_elsayed
    .venv\\Scripts\\python -m tests.generar_fixtures --solo B_base

OJO: los casos de ciclo corren resolver() EN FRIO (caches vacios) para medir lo
mismo que un worker real; A_elsayed tarda ~40 s y B_base ~4 min.
"""
import argparse
import json
import os
import platform
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import kalina as k
import nh3h2o as ng

AQUI = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(AQUI, "fixtures")

# Mismos casos que benchmark/bench.py (PAR_A y BASE_B), copiados para no
# acoplar los tests al modulo de benchmark.
PAR_A = dict(P_alta=1.5, P_baja=0.2735498695811723, w_b=0.55, T_f=373.0,
             T_amb=283.0, eta_t=0.80, eta_p=0.80,
             cierre_HRVG=dict(dT_app=4.0), cierre_reg=dict(dT_pp=4.0),
             cierre_cond=dict(dT_app=4.0))

BASE_B = dict(P_alta=3.0, P_baja=0.4, w_b=0.50, T_f=623.15, T_amb=300.033,
              eta_t=0.85, eta_p=0.85,
              cierre_HRVG=dict(eps=0.85), cierre_reg=dict(eps=0.75),
              cierre_cond=dict(eps=0.80))

# Estados representativos del dominio del proyecto (4-30 bar, 300-630 K,
# w = 0.4-0.7 masica) y puntos de equilibrio biphasico para el motor.
ESTADOS = [
    dict(T=400.0, P=3.0, w=0.50),    # liquido subenfriado
    dict(T=500.0, P=3.0, w=0.50),    # vapor
    dict(T=350.0, P=0.4, w=0.50),    # dentro de la campana (bifasico esperado)
    dict(T=373.0, P=1.5, w=0.55),    # punto de Elsayed et al. (2013)
    dict(T=450.0, P=3.0, w=0.40),    # borde del dominio
    dict(T=320.0, P=0.4, w=0.60),
]

FLASHES = [
    dict(T=340.0, P=0.4),
    dict(T=380.0, P=2.0),
    dict(T=450.0, P=3.0),
]


def meta():
    import numpy, scipy, iapws
    return dict(
        fecha=datetime.now().isoformat(timespec="seconds"),
        python=platform.python_version(),
        plataforma=platform.platform(),
        numpy=numpy.__version__, scipy=scipy.__version__,
        iapws=iapws.__version__)


def guardar(nombre, dato):
    os.makedirs(FIXTURES, exist_ok=True)
    ruta = os.path.join(FIXTURES, nombre)
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(dato, f, indent=1, ensure_ascii=False, allow_nan=True)
    print(f"  {os.path.basename(ruta)}  ({os.path.getsize(ruta):,} B)")


def generar_motor():
    print("MOTOR: estado() y flash_TP()")
    casos = []
    for c in ESTADOS:
        casos.append(dict(**c, salida=k.estado(c["T"], c["P"], c["w"])))
    for c in FLASHES:
        casos.append(dict(**c, flash=k.flash_TP(c["T"], c["P"])))
    guardar("motor.json", dict(meta=meta(), casos=casos))
    return casos


def resolver_caso(nombre, par):
    print(f"CICLO {nombre}: resolver() en frio (caches vacios)...", flush=True)
    est, ind = k.resolver(dict(par))
    est_l = [est[i] for i in range(1, 11)]
    crit = k.criterios(est, ind, par)
    crit_l = [[cid, cumple, sev, txt] for (cid, cumple, sev, txt) in crit]
    etiqueta, violados = k.clasificar(crit)
    guardar(f"ciclo_{nombre}.json", dict(
        meta=meta(), par=par, est=est_l, ind=ind, crit=crit_l,
        etiqueta=etiqueta, violados=violados))
    return est, ind


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo", choices=["motor", "A_elsayed", "B_base"])
    a = ap.parse_args()

    if a.solo in (None, "motor"):
        generar_motor()
    if a.solo in (None, "A_elsayed"):
        resolver_caso("A_elsayed", PAR_A)
    if a.solo in (None, "B_base"):
        resolver_caso("B_base", BASE_B)
    print("Listo. Revisa git diff antes de commitear los fixtures.")


if __name__ == "__main__":
    main()