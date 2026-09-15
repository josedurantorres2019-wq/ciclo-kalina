"""Regresion del ciclo KCS-11 completo: resolver() + criterios + clasificar.

Congela dos casos completos (metodo A de Elsayed 2013 y caso base degenerado
del proyecto). Cada uno tarda 40-250 s en frio: estan marcados como @lento y
se corren explicitamente con `pytest -m lento`.
"""
import pytest

import kalina as k
from tests.helpers import cargar, mismos


def _checkear(nombre):
    fx = cargar(f"ciclo_{nombre}.json")
    est, ind = k.resolver(dict(fx["par"]))
    mismos([est[i] for i in range(1, 11)], fx["est"], camino=f"[{nombre}] est")
    mismos(ind, fx["ind"], camino=f"[{nombre}] ind")
    crit = k.criterios(est, ind, fx["par"])
    mismos([[cid, cumple, sev, txt] for (cid, cumple, sev, txt) in crit],
           fx["crit"], camino=f"[{nombre}] crit")
    etiqueta, violados = k.clasificar(crit)
    assert etiqueta == fx["etiqueta"], f"[{nombre}] etiqueta: {etiqueta}"
    assert violados == fx["violados"], f"[{nombre}] violados: {violados}"


@pytest.mark.lento
def test_ciclo_A_elsayed():
    _checkear("A_elsayed")


@pytest.mark.lento
def test_ciclo_B_base():
    _checkear("B_base")