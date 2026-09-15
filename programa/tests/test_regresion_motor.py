"""Regresion del motor de propiedades: estado() y flash_TP() congelados.

Protege la capa mas baja: si tocamos rho_TPx / P_of / flash_TP y cambia un
punto (T, P, w) del dominio, este test se pone rojo. Rapido (< 10 s).
"""
import pytest

import kalina as k
from tests.helpers import cargar, mismos

FX = cargar("motor.json")


@pytest.mark.parametrize("idx",
                         [i for i, c in enumerate(FX["casos"]) if "salida" in c],
                         ids=lambda i: f"estado{i:02d}")
def test_estado(idx):
    c = FX["casos"][idx]
    actual = k.estado(c["T"], c["P"], c["w"])
    mismos(c["salida"], actual, camino=f"estado({c['T']},{c['P']},{c['w']})")


@pytest.mark.parametrize("idx",
                         [i for i, c in enumerate(FX["casos"]) if "flash" in c],
                         ids=lambda i: f"flash{i:02d}")
def test_flash(idx):
    c = FX["casos"][idx]
    actual = k.flash_TP(c["T"], c["P"])
    mismos(c["flash"], actual, camino=f"flash_TP({c['T']},{c['P']})")
    if actual is not None:
        assert 0 < actual[0] < actual[1] < 1


@pytest.mark.parametrize("idx",
                         [i for i, c in enumerate(FX["casos"]) if "salida" in c],
                         ids=lambda i: f"invariante{i:02d}")
def test_estado_invariantes(idx):
    """Invariantes termodinamicos basicos que deben cumplirse SIEMPRE, no
    solo frente al fixture: h y s crecen con T a P y w fijos (monotonia que
    asumen estado_de y los lazos)."""
    c = FX["casos"][idx]
    e = k.estado(c["T"], c["P"], c["w"])
    e_mas = k.estado(c["T"] + 5.0, c["P"], c["w"])
    assert e_mas["h"] > e["h"], f"h no crece con T en {c}"
    assert e_mas["s"] > e["s"], f"s no crece con T en {c}"