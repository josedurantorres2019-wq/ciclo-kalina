# Plan de optimizacion de rendimiento — ciclo Kalina

Objetivo: reducir el tiempo de un caso (~40 s en frio hoy) manteniendo los
mismos resultados dentro de tolerancia numerica. La red de seguridad ya existe:
`tests/` (golden master sin mocks). Todo cambio pasa por ella.

## Baseline medido (2026-09-15, maquina del usuario) — ANTES de las fases

- A_elsayed en frio: **35.3 s** (neighbor 34.7 s, warm 7.0 s) [pared 88.1 s]
- Distribucion (cProfile): `MEoS._phir` 73.0 % propio (32.3 s); `MEoS.__init__`
  (construir IAPWS95/NH3) 7.0 % acum. — **159 707 construcciones por caso**;
  `P_of` 88.1 %; `rho_TPx` 89.1 %; `perfil_HRVG` 25.2 % (1 llamada, 11.1 s).
- 159 974 llamadas a `_phir` por caso.
- `pytest -m "not lento"` 5.4 s (15 tests) · `pytest -m lento` 253.0 s (2 tests).

## Resultado (fases 1-2 implementadas, 2026-09-15)

Mismas formulas e inversiones; los 17 tests verdes y `validacion_g4_01.py`
conserva `max_dev = 3.85e-03 %` (fisica intacta).

| Metrica | Antes | Despues | Delta |
|---|---:|---:|---:|
| A_elsayed FRIO (bench pared) | 35.3 s | **23.3 s** | **-34 %** |
| A_elsayed caliente | 7.0 s | 5.4 s | -23 % |
| A_elsayed vecino (h1_semilla) | 34.7 s | 21.2 s | -39 % |
| iapws MEoS.__init__ (puros) por caso | 159 707 | **575** | -99.6 % |
| llamadas _phir (puros) por caso | 159 974 | 117 962 | -26 % |
| perfil_HRVG dentro de resolver() | 1 (25 % del perfil) | 0 | fuera del camino caliente |
| `pytest -m "not lento"` | 5.4 s | 5.2 s | -4 % |
| `pytest -m lento` (2 ciclos frios) | 253.0 s | 246.0 s | -2.8 % |

Reparto post-implementacion (cProfile, A_elsayed frio): `MEoS._phir` 77.2 %
propio; `P_of` 90.7 %; `rho_TPx` 91.6 %. El cuello del motor sigue siendo el
cuerpo de `_phir` (fase 5a).

Atribucion honesta de la ganancia: la Fase 2 (C_g,min fuera de resolver) es el
grueso (~11 de los 12 s del frio): el perfil del HRVG era el 2do costo del caso
y ademas condicionaba la robustez de los fixtures a cada cambio del perfil. La
Fase 1 (instancias constantes) quita las ~159 000 construcciones de MEoS por
caso, que en el perfil eran ~7 % (0.30 s propio + 2.8 s acumulado): NO el 73 %
(ese % es el CUERPO de `_phir`, al que las instancias constantes no afectan —
la estimacion original del plan confundio ambos).

Nota de reproducibilidad: el caso B_base (degenerado, cerca del locus critico
donde el flash se estanca) NO es bit-reproducible entre corridas (~1e-10..1e-11
relativo). Se comprobo que regenerando B_base con el codigo ORIGINAL tampoco se
reproduce su fixture previo: el ruido es run-to-run preexistente de la maquina,
no un efecto de estas fases. Dos corridas del MISMO codigo nuevo dan el mismo
fixture (solo cambia la fecha de la metadata). Toda la red usa
rtol=atol=1e-7, varias ordenes por encima del ruido.

## Estrategia

Fases independientes: cada una se fusiona por separado, con medicion antes y
despues. Ordenadas de menor a mayor riesgo/esfuerzo.

## Fase 1 — 1b: construir `IAPWS95()` y `NH3()` una sola vez (73 % del tiempo) — HECHA

> Nota post-implementacion: el 73 % del perfil es el CUERPO de `_phir`, no la
> construccion de los puros. Lo que esta fase elimina son las ~159 000
> construcciones por caso (~7 % del perfil acumulado; 0.30 s propio). La
> medicion real quedo en la seccion Resultado.

En `nh3h2o.py`, el parche de `_phir` (lineas 33-63) recompila el fuente de
`H2ONH3._phir`. El fuente original (iapws 1.5.5) construye en cada llamada:

    water = IAPWS95()        ; lineas 68-69 del fuente
    ammonia = NH3()          ; lineas 71-72 del fuente

Esas dos instancias SON CONSTANTES: solo se usan para `_phir(tau, delta)`,
que depende solo de (tau, delta), no del estado interno. Plan:

1. En `_aplicar_correcciones()`, crear dos globales `_WATER = IAPWS95()` y
   `_AMMONIA = NH3()` y agregarlas al dict `ns`.
2. Agregar a `_SUBS["_phir"]` los reemplazos exactos:
   - `"        water = IAPWS95()\n"` -> `"        water = _WATER\n"`
   - `"        ammonia = NH3()\n"` -> `"        ammonia = _AMMONIA\n"`
3. Verificar que `src.replace` de esas cadenas ocurra 1 vez (falla la
   guarda si el paquete instalado difiere).

Esperado: ~70 % menos de tiempo total por caso. Debe reproducir los numeros
exactos (mismas llamadas, mismas formulas). `validacion_g4_01.py` y los tests
deben pasar SIN regenerar fixtures.

## Fase 2 — 3b: `C_g,min` bajo demanda — HECHA

Hoy `resolver()` (kalina.py lineas 718-727), cuando la fuente es isoterma
(`C_g` sin declarar, el caso tipico), corre SIEMPRE `perfil_HRVG` +
`C_g_minimo_HRVG` para rellenar `ind['Cg_min_HRVG']`. Ese dato es solo
INFORMATIVO: solo lo leen `criterios()` (linea 799) y `supuestos()` (linea
956) para armar un texto. No condiciona la convergencia.

Plan:
1. En `resolver()`, rama isoterma: sacar el calculo; dejar
   `ind['Cg_min_HRVG'] = None`.
2. En `supuestos()`/`criterios()`, cuando necesiten el texto, calcularlo bajo
   demanda (helper `_cg_min_hrvg(est, ind, par)` que llama `perfil_HRVG` +
   `C_g_minimo_HRVG` con los datos ya disponibles en `est`/`ind`/`par`).
3. Como los textos se generan igual, `criterios()` y `supuestos()` dan los
   mismos resultados; PERO `ind['Cg_min_HRVG']` cambia de 15.18/8.23 a `null`
   en los fixtures de ciclo. Ese cambio es INTENCIONAL:
   - regenrar `ciclo_A_elsayed.json` y `ciclo_B_base.json` SOLO para reflejarlo,
   - `git diff` del fixture debe mostrar exactamente ese cambio y nada mas,
   - y además `ciclo_*.json` queda mas liviano (el json pasa de ~6 KB a ~1 KB).
4. El texto "C_g,min = X" debe seguir apareciendo identico en supuestos (la
   prueba: correr `supuestos()` por consola antes/despues y comparar).

Esperado: ahorra la primera evaluacion del perfil del HRVG (~1-2 s? medir).
Ademas hace que los fixtures de ciclo no regen por cada cambio en el perfil.

## Fase 3 — 1c: Newton salvaguardado en `rho_TPx`

`nh3h2o.py` lineas 89-111: el bubbleo de densidad (raiz de P(rho)=P dada, con
arbitrario de fase) usa `_bracket_log` + `brentq`, ~17 evaluaciones de P_of
por raiz. Con la derivada ∂P/∂ρ disponible (analitica o en diferencias sobre
la misma P_of), un Newton salvaguardado con reinicio por bracket llega en
~4-6 evaluaciones.

Requisito: el Newton DEBE dar el mismo rho final (las mismas raices del
fisico). Con la red: si `motor.json` o los tests de ciclo se ponen rojos por
ruido numerico (< 1e-5), es aceptable regenerar SOLO tras validar que
`validacion_g4_01.py` da el mismo max_dev. Documentar la desviacion maxima
observada entre brentq y Newton en un barrido de (T, P, w).

## Fase 4 — 3a: Newton con c_p en `estado_de`

`kalina.py` lineas 161-203: invierte h(T) → T con brentq y bracket expandido
(~12 evaluaciones de `estado()`). Reemplazarlo por Newton salvaguardado que
arranca con la derivada c_p = ∂h/∂T estimada por `estado(T+δ)` (una sola
evaluacion extra por iteracion). Mismo criterio que Fase 3.

## Fase 5 — (fuera de alcance por ahora) 1a/2a

- 1a: `P_of` especializado que calcule solo fir/fird/firdd/firx (lo que usan
  los solvers) sin firt/firtt/firdt. Ganancia ~4x adicional sobre P_of.
- 2a: Jacobiano analitico en el flash_TP para resolver el par (xL, xV).

Se delegan despues si las fases 1-4 no llegan al objetivo (caso ~5 s).

## Protocolo de verificacion obligatorio (en cada fase)

1. `pytest -m "not lento"` (motor, 7 s) — debe quedar verde.
2. `pytest -m lento` (ciclos, 5 min) — verde; si rompe por ruido < 1e-5,
   regenerar fixtures SOLO si `validacion_g4_01.py` conserva max_dev.
3. `python validacion_g4_01.py` — max_dev igual (3.85e-03 %) salvo ruido.
4. `python bench.py --escenarios A_elsayed` frio — medir nuevo tiempo.
5. Perfil con cProfile (`python -m cProfile`) sobre A_elsayed y reportar el
   nuevo % de `_phir`, `P_of`, `rho_TPx`.
6. Tabla de tiempos por fase en este documento.

## Reglas

- NO tocar la fisica: mismas formulas, mismas inversiones, mismas guardas.
- NO reescribir fixtures en bloque: cada regeneracion se justifica con la
  causa y se revisa con `git diff`.
- NO subir archivos binarios ni `.venv`.
- Commits convencionales cortos por fase (si el usuario los pide).