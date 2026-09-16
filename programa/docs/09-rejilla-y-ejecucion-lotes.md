# 09 — Rejilla de casos y ejecución por lotes

**Archivos:** [`motor_ui.py`](../motor_ui.py), [`worker_lote.py`](../worker_lote.py), [`worker_caso.py`](../worker_caso.py), [`app.py`](../app.py)

## Qué resuelve

Convierte la configuración de la interfaz (cada variable FIJO o BARRIDO) en una
lista de casos, los resuelve de forma robusta ante casos que se cuelgan, los
presenta y los registra en Excel. Además implementa los modos SELECCION y
FRONTERA sobre los resultados.

## Fundamento

### Continuación numérica (arranque tibio)

Si dos casos vecinos de la rejilla difieren poco en sus parámetros, su solución
h₁ también difiere poco (la solución es continua fuera de los cambios de
régimen). Usar h₁ del caso anterior como semilla del siguiente reduce el
intervalo de búsqueda de Brent:

$$
h_1^{(k+1)} \in \left[h_1^{(k)} - \Delta,\; h_1^{(k)} + \Delta\right],\qquad \Delta = 2\cdot 3^{j}\ \text{kJ/kg}
$$

Lo mismo ocurre con las tablas internas (`_flash`, `_ultimaT`): se "calientan" y
aceleran los casos siguientes **dentro del mismo proceso**.

### Por qué procesos y no hilos

- Los cachés de `kalina.py` son globales de módulo, sin locks.
- Una llamada colgada dentro de código numérico no se puede cancelar desde un
  hilo; matar el proceso sí.
- `subprocess` evita que Windows (método *spawn* de `multiprocessing`) reimporte
  la app de Streamlit en el hijo.

## Algoritmos

### A. Rejilla — [`construir_grid`](../motor_ui.py#L111)

Producto cartesiano de los valores de cada variable
(`np.linspace(min, max, n)` si es BARRIDO). N casos = ∏ nᵢ.

### B. Ejecución — [`resolver_grid`](../motor_ui.py#L193)

```
pendientes = todos
repartir pendientes en K bloques contiguos (K = n_workers)
lanzar un worker_lote.py por bloque (la rejilla en paralelo)
por cada worker:
    hilo lector stdout → eventos a una cola compartida
    para cada caso esperado de ese worker:
        esperar linea con timeout (por worker, el deadline mas proximo de todos)
        ├─ llegó → registrar, avanzar al siguiente caso del bloque
        ├─ timeout → matar SOLO ese worker, marcar TIEMPO_AGOTADO,
        │            re-marcar el resto de su bloque como pendiente
        └─ fin inesperado → marcar ERROR con stderr, re-marcar el resto
siempre que haya un worker libre y casos pendientes, lanzar el siguiente
```

El worker ([`worker_lote.main`](../worker_lote.py#L57)) resuelve su bloque en orden,
encadena `h1_semilla` entre casos OK y escribe una línea JSON con `flush` por caso.

La paralelización (`n_workers`, default auto = núcleos − 1) solo pierde el caché de
un worker cuando un caso CONCRETO de ese bloque se cuelga: los demás siguen
intactos. La continuidad numérica se preserva porque la rejilla se parte en
bloques contiguos (la semilla h₁ encadena dentro de cada bloque). Rejillas de
< 8 casos corren serial (1 worker): el arranque frío pesa más que el paralelismo
en rejillas chicas.

### C. SELECCION — [`aplicar_seleccion`](../motor_ui.py#L361)

Filtra VALIDO / VALIDO con advertencia, excluye O1, ordena por W_net y η
descendente. Solo se permite con BARRIDO en P_alta, P_baja, x_b, m_b.

### D. FRONTERA — [`calcular_frontera`](../motor_ui.py#L374)

Agrupa por las demás variables barridas, ordena por el parámetro de cierre y
reporta cada par de puntos consecutivos donde cambia "es Kalina". No interpola:
el separador es discontinuo (paso de bifásico a monofásico en el estado 2), así
que una interpolación daría una frontera falsa.

### E. Registro — [`guardar_en_excel`](../motor_ui.py#L343)

Lee el libro completo, concatena las filas nuevas y lo reescribe.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Paralelismo con K workers** | ✅ **HECHO**. Cada bloque contiguo corre en su propio `worker_lote.py` (`resolver_grid(n_workers=...)`); timeout y reinicio por worker, stderr drenado en thread (sin bloqueo de pipe) | Hecho |
| **Orden serpenteante** | `itertools.product` salta del último valor de una fila al primero de la siguiente: la semilla h₁ queda lejos justo ahí. Recorrer en *boustrophedon* (invertir el sentido de la variable interna en cada fila) mantiene todos los vecinos adyacentes | Bajo |
| **Semilla del vecino, no del anterior** | En rejillas 2D+ el caso previo en la lista no siempre es el vecino más próximo. Guardar h₁ por índice de rejilla y sembrar con el vecino ya resuelto más cercano | Medio |
| **Riesgo de bloqueo por `stderr`** | ✅ **HECHO**: `resolver_grid` drena `stderr` en un thread por worker (buf en memoria, no `PIPE` sin leer), así un exceso de mensajes no bloquea al proceso ni produce TIEMPO_AGOTADO falso | Hecho |
| **Excel incremental** | `guardar_en_excel` es O(filas totales) por guardado y crece sin fin. Usar `openpyxl` en modo *append* sobre la hoja, o un CSV/Parquet de trabajo que se exporte a Excel bajo demanda | Bajo |
| **Relanzar tras crash** | Tras un fin inesperado se marca ERROR solo el caso en curso y se relanza; si el crash es sistemático (p. ej. importación fallida) se relanza N veces. Abortar si el worker muere sin producir ninguna línea dos veces seguidas | Bajo |
| **`pendientes.remove(i)`** | O(n) por caso, O(n²) total; irrelevante hasta miles de casos, pero un `set` o un índice lo resuelve | Trivial |
| **Refinamiento de FRONTERA** | La frontera se localiza con la resolución de la rejilla. Una segunda pasada automática que biseque solo los tramos con cruce (sin interpolar, resolviendo casos reales) da la frontera a la precisión deseada con log₂ de las evaluaciones | Medio |
| **`worker_caso.py`** | Ya no lo usa la interfaz (`resolver_caso` no se llama desde `app.py`). Mantenerlo solo si se quiere el aislamiento total por caso; si no, eliminarlo evita duplicar lógica con `worker_lote.py` | Trivial |
