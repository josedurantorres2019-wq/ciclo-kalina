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

## Resultado (fases 3-4 implementadas, 2026-09-15)

Baseline de esta sesion (post fases 1-2, misma maquina): A_elsayed FRIO 22.4 s;
`pytest -m "not lento"` 7.16 s; `pytest -m lento` 238.3 s (2 verdes); P_of por
caso 55 283; brentq por caso 3 354; `validacion_g4_01.py` byte-identico.

| Metrica | Fases 1-2 | +Fase 3 | +Fases 3-4 |
|---|---:|---:|---:|
| A_elsayed FRIO (bench pared) | 22.4 s | 23.5 s | **18.8 s (-16 %)** |
| A_elsayed vecino (h1_semilla) | 20.3 s | 19.9 s | **16.8 s (-17 %)** |
| A_elsayed caliente | 5.3 s | 6.5 s | 5.4 s |
| pared total (frio+diag+caliente+vecino) | 66.1 s | 68.6 s | **60.9 s (-8 %)** |
| llamadas brentq (scipy) por caso | 3 354 | 39 | **4** |
| llamadas P_of por caso | 55 283 | 57 599 | **45 396 (-18 %)** |
| llamadas rho_TPx por caso | 3 316 | 3 503 | **2 734** |
| llamadas estado() por caso | 386 | 386 | **316** |
| llamadas MEoS._phir (puros) por caso | 117 962 | 122 970 | **96 994 (-18 %)** |
| `pytest -m "not lento"` | 7.16 s | 5.17 s | 5.36 s (-25 %) |
| `pytest -m lento` (2 ciclos frios) | 238.3 s | 246.3 s | **215.0 s (-10 %)** |

`validacion_g4_01.py` queda BYTE-IDENTICO al baseline (max_dev 3.85e-03 %).
Fixtures de ciclo regenerados dos veces (fase 3 y fase 4); `git diff` mostro
solo cifras finas (max desviacion relativa 3.6e-09 en A_elsayed, 1.3e-09 en
B_base, frente al criterio 1e-5) y fecha de metadata; `motor.json` NO se
regenero (la red rapida quedo verde en las dos fases). `etiqueta` (VALIDO) y
`violados` sin cambios.

### Lectura honesta de la medicion

- **Fase 3 (Newton en rho_TPx) es casi NEUTRA en pared para el ciclo**: la
  premisa del plan ("~12-17 evaluaciones de brentq por raiz") NO se cumple.
  Sobre estas funciones monotonas brentq hace ~7.5 evals/raiz, y el costo real
  lo domina la busqueda de bracket (expansiones *1.06/*0.97): 9.2 evals/raiz
  de media. El Newton salvaguardado gasta ~7.3 evals/raiz en el apriete
  (2 evals por iteracion + evaluacion del extremo): empate tecnico con brentq.
  Su beneficio CONCRETO: saca a scipy brentq del hot path (3 354 -> 4 llamadas)
  y acelera el suite rapido del motor (7.2 -> 5.2 s, -25 %).
- **Fase 4 (Newton en estado_de) es LA ganancia**: estado_de baja de ~11-12
  `estado()` por inversion a ~5 (2-3 iteraciones de Newton con c_p efectiva).
  Eso recorta 316 vs 386 `estado()`, 2 734 vs 3 316 raices rho_TPx y el P_of
  total un 18 %. Traduccion en pared: FRIO -16 %, vecino -17 %, lento -10 %.
- Efecto colateral medido: la perturbacion de convergencia (~1e-10) cambia las
  rutas de fsolve del flash y, en la medicion SOLO fase 3, subio P_of total
  (+4 %) frente al baseline; con fase 4 el neto quedo globalmente -18 %.

Nota de reproducibilidad: el caso B_base (degenerado, cerca del locus critico
donde el flash se estanca) NO es bit-reproducible entre corridas (~1e-10..1e-11
relativo). Se comprobo que regenerando B_base con el codigo ORIGINAL tampoco se
reproduce su fixture previo: el ruido es run-to-run preexistente de la maquina,
no un efecto de estas fases. Dos corridas del MISMO codigo nuevo dan el mismo
fixture (solo cambia la fecha de la metadata). Toda la red usa
rtol=atol=1e-7, varias ordenes por encima del ruido.

## Resultado (fase 1a implementada, 2026-09-16)

`P_of` especializado que evalua SOLO `fird = dPhi_r/ddelta` de los puros y del
termino de mezcla (en vez de las 6-7 derivadas de `iapws _phir`). Mismo orden
de operaciones que el original: verificado BIT-IDENTICO en 20 000 barridos de
(rho, T, x). La fisica no se toca: `validacion_g4_01.py` da el mismo
`max_dev = 3.85e-03 %` y `pytest -m lento` pasa SIN regenerar fixtures.

| Metrica | Fases 3-4 | +Fase 1a | Delta |
|---|---:|---:|---:|
| A_elsayed FRIO (bench pared) | 18.8 s | **4.8 s** | **-75 %** |
| A_elsayed caliente | 5.4 s | **1.2 s** | -78 % |
| A_elsayed vecino (h1_semilla) | 16.8 s | **4.1 s** | **-76 %** |
| B_base FRIO | ~107 s (2015 lento/2) | **46.8 s** | **-80 %** |
| B_base vecino | ~87 s | **34.3 s** | -61 % |
| diagnostico (criterios+supuestos) | 19.1 s | **3.8 s** | -80 % |
| `pytest -m "not lento"` | 5.36 s | **2.19 s** | -59 % |
| `pytest -m lento` (2 ciclos frios) | 215.0 s | **59.0 s** | **-73 %** |
| llamadas MEoS._phir (puros) por caso | 96 994 | **6 202** | -94 % |
| % del perfil en MEoS._phir (propio) | 76.1 % | **20.7 %** | |
| perfilado total (A_elsayed) | 26.2 s | **6.3 s** | -76 % |

Reparto post-implementacion (cProfile, A_elsayed frio): `fird` puro propio
42.9 %; `MEoS._phir` (ya solo desde `prop()`) 20.7 %; `_fird_departure` 8.1 %;
`math.exp` 6.9 %. El antiguo monstruo (70+ %) quedo reducido a 1/4.

Cambios (todos en `nh3h2o.py`): `_mk_fird_puro()` genera el `fird` podado de
cada puro desde sus `_constants` (IAPWS-95 con polinomios+exp+gaussiano+no
analitico; NH3 solo polinomio+exp); `_fird_departure()` replica el `fird` del
termino de mezcla de `_Dphir` (misma tabla de coeficientes, mismo orden);
`_delta_fird()` combina (1-x)·fird_agua + x·fird_amon + fird_mezcla; `P_of`
lo consume. Se importa `math.exp` a nivel de modulo (bits identicos a
`numpy.exp`, medido).

Nota de atribucion honesta: con la capa de propiedades ~4x mas barata, el
reparto se achato. El otro hot path es ahora la SECUENCIA de llamadas
(flash→F→rho_TPx→P_of): sus cuentas por caso NO cambiaron (45 396 P_of,
2 734 prop), solo el costo unitario. Para bajar MAS por caso hay que bajar el
NUMERO de evaluaciones: Jacobiano explicito en flash_TP (fase 2a) para
fsolve -> menos F por flash, y semilla interpolada (2b). Para un barrido de
muchos ciclos, la palanca ortogonal es paralelizar la rejilla en K procesos
worker_lote (doc 09, HECHO 2026-09-16: `resolver_grid(n_workers=...)`).

Dato de barrido: B_base vecino sigue costando 34.3 s (solo -61 % vs 2015):
en un barrido con T_f creciente casi ningun flash acierta el cache (que
guarda por T exacta) y cada caso paga la inversion completa. La mejora de
vecino es promesa de 2a/2b, no de la capa de propiedades.

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

## Fase 3 — 1c: Newton salvaguardado en `rho_TPx` — HECHA

> Nota post-implementacion: la ganancia esperada NO aparece en pared (el ciclo
> quedo neutro; ver tablas de arriba) porque brentq ya convergia en ~7.5
> evals/raiz sobre estas funciones monotonas y el bracket domina (~9.2 evals).
> El Newton salvaguardado elimina igual el brentq del camino caliente
> (3 354 -> 4 llamadas) y acelera el suite rapido del motor.

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

## Fase 4 — 3a: Newton con c_p en `estado_de` — HECHA

> Nota post-implementacion: es la fase que entrega el grueso de la ganancia.
> `estado_de` pasa de ~11-12 `estado()` por inversion a ~5; recorta FRIO de
> 22.4 a 18.8 s (-16 %), vecino de 20.3 a 16.8 s (-17 %) y el lento de 238 a
> 215 s (-10 %). Implementacion: `_newton_salva_T` (kalina.py), mismo contrato
> de salvaguardas que la Fase 3; x0 = ultima T de `_ultimaT` si cae en el
> bracket, si no el extremo de menor residuo; convergencia |ΔT| <= TOL_T.

`kalina.py` lineas 161-203: invierte h(T) → T con brentq y bracket expandido
(~12 evaluaciones de `estado()`). Reemplazarlo por Newton salvaguardado que
arranca con la derivada c_p = ∂h/∂T estimada por `estado(T+δ)` (una sola
evaluacion extra por iteracion). Mismo criterio que Fase 3.

## Fase 1a — `P_of` especializado (solo fird) — HECHA

> Implementada 2026-09-16. `P_of` ahora evalua solo `dPhi_r/ddelta` de los
> puros (`_mk_fird_puro`, podando MEoS._phir) y del termino de mezcla
> (`_fird_departure`, podando _Dphir), con el mismo orden de operaciones
> (BIT-IDENTICO en 20 000 barridos). FRIO 18.8 -> 4.8 s (A), 107 -> 46.8 s
> (B, -80 %); MEoS._phir baja de 96 994 a 6 202 llamadas/caso. Resultado
> completo en la seccion "Resultado (fase 1a implementada, 2026-09-16)".

## Fase 2b — Semilla interpolada en `flash_TP` — HECHA

> Implementada 2026-09-16 (kalina.py, `flash_TP`). Cuando la tabla `_flash[P]`
> ya tiene soluciones a ambos lados de T, interpola x_L/x_V linealmente entre
> los dos vecinos en vez de usar el mas proximo (la campana zeotropica es casi
> lineal en T a P fija). Efecto medido (bench.py, 20260916_114435):
>
> | Escenario | post-1a | +2b | Ganancia |
> |---|---|---:|---:|
> | A_elsayed FRIO | 4.8 s | 4.75 s | -1 % |
> | A_elsayed vecino | 4.1 s | 3.9 s | -5 % |
> | B_base FRIO | 46.8 s | 37.7 s | **-19 %** |
> | B_base vecino | 34.3 s | 28.1 s | **-18 %** |
>
> Las evaluaciones de F por flash bajan de ~4.7-7.6 a ~3.7 (A_elsayed). En
> B_base se mantienen ~296 reintentos con Raoult (fsolve 1817 vs 1521 flashes).
> La ganancia es mayor donde la tabla es grande (barridos).

## Fase 5 — (fuera de alcance por ahora) 2a

- 2a: Jacobiano analitico en `flash_TP` (Newton 2x2 explicito o Jacobiano
  via derivacion implicita de rho_TPx) para resolver el par (xL, xV).
  fsolve hoy hace ~3.7 evaluaciones de F por flash (2b implementada); un
  Newton con Jacobiano explicito podria bajar a ~2. Es la siguiente palanca
  por caso, con esfuerzo ALTO.

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