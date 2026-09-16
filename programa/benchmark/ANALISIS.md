# Análisis del benchmark — qué optimizar primero

Corrida: 2026-09-14, Python 3.13.7, Windows 11.
Datos completos: [resultados/20260914_122556/reporte.md](resultados/20260914_122556/reporte.md)
(JSON y `.prof` en la misma carpeta; `python -m pstats <archivo>.prof` para explorarlos).

## 1. Resultado en una frase

**El 92 % del tiempo de cualquier caso se va en evaluar la presión P(ρ) dentro
de la búsqueda de densidad**, y tres cuartas partes de todo el tiempo están en
una sola función de `iapws`: la energía de Helmholtz residual de los componentes
puros (`MEoS._phir`), que calcula todas las derivadas aunque `P_of` solo necesita
una. El patrón es idéntico en los 5 escenarios (±1 %). Lo que cambia entre
escenarios es **cuántas veces** se llega ahí, no **dónde** se gasta.

## 2. Tiempos de pared (sin perfilador)

| Escenario | Qué ejercita | Frío [s] | Caliente [s] | Vecino [s] |
|---|---|---:|---:|---:|
| A_elsayed | Método A, sin lazo exterior | 39.5 | 7.3 | 34.7 |
| B_kalina | Lazo exterior + lazo frío, régimen Kalina | 200.3 | 47.5 | 162.9 |
| B_base | Caso base del proyecto (degenerado) | 240.8 | 43.6 | 174.9 |
| B_base_Cg | + limitador de pinzamiento del HRVG | 291.1 | 62.9 | 200.0 |
| B_base_Ccold | + perfil del condensador | **falla** (§6) | | |

- **Vecino** (T_f + 2 K con `h1_semilla`, la situación real de un barrido) cuesta
  el **72–88 % del frío**. El arranque tibio casi no ayuda: un barrido de N puntos
  cuesta ≈ N × frío.
- **Caliente** (exactamente el mismo caso otra vez) sigue costando 44–63 s en los
  casos B: los cachés guardan composiciones de equilibrio, pero las densidades y
  las propiedades se recalculan siempre.
- **Memoria:** 76 MB de pico en todos los casos; los cachés llegan a ~2 000
  soluciones de flash. **La memoria no es un problema.**
- **Sobrecosto del perfilador:** +12 % a +24 % sobre el tiempo de pared. Los
  porcentajes del perfil son fiables para ordenar.

## 3. Arranque

| Paso | Tiempo |
|---|---:|
| Relanzar un worker (`python` + `import kalina`) | 1.23 s |
| … de lo cual `import scipy.optimize` | 0.86 s |
| … parches de `iapws` en `nh3h2o` | 0.005 s |
| Abrir la interfaz (`motor_ui` + `streamlit` + `pandas`) | 1.86 s |
| Primera `estado()` en frío | 0.47 s |
| `estado()` en el mismo punto (flash cacheado) | 0.016 s |
| `estado()` a 1 K de distancia | 0.25 s |
| Primera inversión `estado_de()` | 3.09 s |
| Inversión vecina `estado_de()` | 1.82 s |

**El arranque no es donde está el problema**: 1–2 s frente a 40–300 s por caso.
Lo que sí importa es que **una sola inversión h→T cueste 2–3 s**: un caso B hace
unas 140.

## 4. Dónde se va el tiempo — por capas

% del tiempo total de `resolver()` en frío (rango en los 5 escenarios).

```
resolver                                   100 %
└─ tramo(h1)  (lazo exterior)               65–100 %
   ├─ lazo_frio                             29–84 %
   └─ estado_de  (inversión h/s→T)          60–95 %
      └─ estado                             100 %
         └─ flash_TP → fsolve → F           81–92 %
            └─ rho_TPx  (raíz de densidad)  92 %
               └─ P_of  (P(ρ))              91–92 %
                  └─ H2ONH3._phir (mezcla)  96–97 %
                     ├─ MEoS._phir (puros)  75 % PROPIO   ◄── punto caliente
                     ├─ _Dphir parcheado    14 % PROPIO
                     └─ construir IAPWS95()/NH3()  3–4 %
```

### Números clave que explican el costo

| Métrica | Valor | De dónde sale |
|---|---:|---|
| Evaluaciones de P(ρ) por raíz de densidad | **≈ 17** | `P_of` / `rho_TPx` = 478 660 / 28 351 (B_base) |
| Objetos `IAPWS95()`/`NH3()` construidos | ≈ 1 000 000 por caso | 2 por cada llamada a `_phir` de la mezcla |
| Evaluaciones de F por `fsolve` | 4.7–7.6 | incluye las del Jacobiano por diferencias finitas |
| Llamadas a `fsolve` / llamadas a `flash_TP` | 0.98–1.34 | **el caché de flash casi nunca acierta** (solo con T idéntica) y a veces falla la semilla tibia y se reintenta con Raoult |
| `estado()` por inversión `estado_de` | ≈ 12 | B_base: (1 731 − 40) / 138 |
| Costo del diagnóstico `C_g,min` al final de `resolver` | 11–13 s | `perfil_HRVG` (1 llamada): 27 % en A, 4–5 % en B |
| `criterios()` en método A | 9.3 s | `T_sat` → `bubbleT` completo (criterio F1) |

## 5. Lista priorizada

Techo teórico (ley de Amdahl): si P(ρ) fuera gratis, un caso sería como máximo
1/(1 − 0.92) ≈ **12 veces** más rápido. Las mejoras de capas distintas se
**multiplican** (menos inversiones × menos flashes × menos raíces × P(ρ) más barata).

### Prioridad 1 — capa de propiedades (`nh3h2o.py`, doc 01) · 92 % del tiempo

| # | Cambio | Qué ataca | Esfuerzo | Riesgo |
|---|---|---|---|---|
| 1a | `P_of` especializado que calcule solo Φʳ y ∂Φʳ/∂δ de los puros y del término de mezcla, sin pasar por `MEoS._phir` completo | 75 % + 14 % propio | Alto | Medio: validar con `validacion_g4_01.py` |
| 1b | Construir `IAPWS95()` y `NH3()` **una vez** (hoy 2 por llamada) | 3–4 % | Muy bajo | Bajo |
| 1c | Newton salvaguardado en `rho_TPx` con ∂P/∂ρ | ≈ 17 → ≈ 4–5 evaluaciones por raíz | Medio | Bajo |
| 1d | Semilla de densidad desde la última raíz de (T, P, x, fase) cercana | Acotamiento geométrico | Bajo | Bajo |

1b es lo más barato de hacer; 1c es probablemente la mejor relación
esfuerzo/ganancia; 1a es la de mayor ganancia.

### Prioridad 2 — flash y estado (`kalina.py`, doc 02)

| # | Cambio | Qué ataca | Esfuerzo |
|---|---|---|---|
| 2a | Jacobiano analítico en `flash_TP` (o Newton 2×2 explícito) | ~3.7 evaluaciones de F por flash (2b ya implementada), la mayoría para el Jacobiano | Alto |
| 2b | Semilla **interpolada** entre los dos flashes vecinos de la tabla | Reintentos con Raoult (fsolve > flash en B) | Bajo · **HECHA (2b)**: B_base -19 % FRIO, -18 % vecino |
| 2c | Devolver desde el flash las densidades de la última evaluación y reutilizarlas en `_mono` | `_mono`: 7–15 % | Medio |

### Prioridad 3 — inversión y lazos (`kalina.py`, docs 03–06)

| # | Cambio | Qué ataca | Esfuerzo |
|---|---|---|---|
| 3a | Newton/secante con c_p en `estado_de` | ≈ 12 `estado()` por inversión; `estado_de` es 60–95 % | Medio |
| 3b | **`C_g,min` bajo demanda**, como ya se hace con `C_cold_min` | 11–13 s fijos por `resolver()` sin C_g (27 % en A) | Muy bajo |
| 3c | Semilla de h₁₀ entre llamadas + Aitken en `lazo_frio` | lazo frío: 29–84 % | Bajo |
| 3d | Perfil del condensador tabulado una vez por (P_baja, w_b) | 39 % con C_cold | Bajo–medio |
| 3e | Perfil del HRVG tabulado una vez por (P_alta, w_b) | 19 % con C_g | Bajo–medio |
| 3f | F1: evitar `bubbleT` en `criterios()` | 9.3 s en método A | Bajo |

### No prioritario

- **Arranque** (1–2 s), **memoria** (76 MB), **Excel/pandas**: irrelevantes
  frente al costo por caso.
- Paralelizar la rejilla (doc 09) multiplica el rendimiento por el número de
  núcleos, pero **no reduce** el costo de un caso: conviene después de las
  prioridades 1–3.

### Orden sugerido

1. **3b y 1b** — minutos de trabajo, ganancia inmediata medible, riesgo casi nulo.
2. **1c** — Newton en densidad.
3. **3a** — Newton en la inversión.
4. **1a / 2a** — derivadas analíticas: la mayor ganancia y el mayor trabajo.

Después de cada paso: `python validacion_g4_01.py` y `python prueba5.py`
(regresión numérica), y `python benchmark/bench.py --escenarios A_elsayed B_base`
para medir la ganancia contra esta línea base.

## 6. Hallazgo colateral: B_base_Ccold falla

`B_base` con `C_cold = 10 kW/K` (sin otros cambios) termina con

```
RuntimeError: equilibrio no resoluble en T=448.3674 K, P=3 MPa, w=0.5000
```

tras ~4 min. El error sale de `estado()` ([kalina.py:114](../kalina.py#L114)) del
**lado de alta presión** (P = 3 MPa, w = w_b), no del condensador: el perfil
del condensador cambia h₉ y con ello los estados que recorre el lazo exterior,
que llega a un punto donde el flash no converge y las funciones de respaldo
(`bubbleP`/`dewP`) dejan P entre burbuja y rocío. Como es `RuntimeError` y no
`NoConverge`, la interfaz lo muestra como **ERROR** y no como **NO_CONVERGE**.
El traceback no quedó guardado en esta corrida (el benchmark ya se corrigió para
guardarlo).
