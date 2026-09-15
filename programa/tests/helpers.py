"""Comparador numerico tolerante y carga de fixtures (golden master).

Los fixtures son la VERDAD congelada del comportamiento actual del motor.
Cualquier cambio en los numeros por encima de la tolerancia es una regresion
(o un cambio intencional que hay que regenerar y revisar con git diff).
"""
import json
import math
import os

AQUI = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(AQUI, "fixtures")


def cargar(nombre):
    with open(os.path.join(FIXTURES, nombre), encoding="utf-8") as f:
        return json.load(f)


def mismos(a, b, rtol=1e-7, atol=1e-7, camino="valor"):
    """Compara a y b recursivamente: numeros con math.isclose (NaN==NaN,
    inf==inf), el resto con igualdad estricta."""
    if isinstance(a, dict) and isinstance(b, dict):
        assert a.keys() == b.keys(), f"{camino}: claves {sorted(a)} vs {sorted(b)}"
        for k in a:
            mismos(a[k], b[k], rtol, atol, f"{camino}.{k}")
        return
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        assert len(a) == len(b), f"{camino}: longitudes {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            mismos(x, y, rtol, atol, f"{camino}[{i}]")
        return
    if isinstance(a, bool) or isinstance(b, bool):
        assert a is b or a == b, f"{camino}: {a!r} != {b!r}"
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if math.isnan(a) and math.isnan(b):
            return
        if math.isinf(a) or math.isinf(b):
            assert a == b, f"{camino}: {a!r} != {b!r}"
            return
        if a == b or math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
            return
        raise AssertionError(
            f"{camino}: {a!r} != {b!r} (rtol={rtol}, atol={atol}, "
            f"diff_rel={abs(a-b)/max(abs(b),atol):.2e})")
    assert a == b, f"{camino}: {a!r} != {b!r}"


def tol(referencia, rtol=1e-7, atol=1e-7):
    """Devuelve tolerancias para el comparador (por ahora fijas; aqui se
    centralizan para ajustarlas de una sola vez si el ruido de maquina las
    exige)."""
    return rtol, atol