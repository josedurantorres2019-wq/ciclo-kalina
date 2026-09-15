# Tests de regresión — golden master

El proyecto no usa mocks: los tests congelan las salidas REALES del motor
(fixtures JSON) y cualquier modificación debe reproducirlas. Es la red de
seguridad antes de optimizar el rendimiento.

## Cómo correr

Desde `programa/` (con el `.venv` activado):

```powershell
# Rapido: motor de propiedades (~10 s)
.venv\Scripts\python -m pytest -v -m "not lento"

# Lento: ciclos completos en frio (~5 min: A_elsayed 40 s + B_base 4 min)
.venv\Scripts\python -m pytest -v -m lento

# Todo
.venv\Scripts\python -m pytest
```

## Regenerar los fixtures

Solo ante un cambio INTENCIONAL de los números (nueva física, parámetros,
versión de iapws). Después de regenerar, `git diff` DEBE mostrar exactamente
el cambio esperado y nada más:

```powershell
.venv\Scripts\python -m tests.generar_fixtures
# o, para casos puntuales:
.venv\Scripts\python -m tests.generar_fixtures --solo motor
.venv\Scripts\python -m tests.generar_fixtures --solo A_elsayed
.venv\Scripts\python -m tests.generar_fixtures --solo B_base
```

## Regla de oro

| Caso | Qué protege | Costo |
|---|---|---|
| `motor.json` | `estado()`, `flash_TP()` en el dominio del proyecto | ~10 s |
| `ciclo_A_elsayed.json` | Método A completo (Elsayed 2013) | ~40 s |
| `ciclo_B_base.json` | Caso base degenerado del proyecto | ~4 min |

Además de estos tests, seguí corriendo la validación externa:
`validacion_g4_01.py` (vs IAPWS) y `prueba5.py` (consistencia A↔B).