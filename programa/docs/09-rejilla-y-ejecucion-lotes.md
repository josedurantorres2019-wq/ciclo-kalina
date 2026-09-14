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

### B. Ejecución — [`resolver_grid`](../motor_ui.py#L192)

```
pendientes = todos
mientras haya pendientes:
    lanzar worker_lote.py con la lista pendiente por stdin (JSON)
    hilo lector: cada línea de stdout → cola
    para cada caso esperado:
        esperar línea con timeout
        ├─ llegó → registrar, reiniciar cronómetro
        ├─ timeout → matar proceso, marcar TIEMPO_AGOTADO, relanzar con el resto
        └─ fin inesperado → marcar ERROR con stderr, relanzar con el resto
```

El worker ([`worker_lote.main`](../worker_lote.py#L57)) resuelve en orden, encadena
`h1_semilla` entre casos OK y escribe una línea JSON con `flush` por caso.

Solo se pierde el caché cuando un caso concreto se cuelga.

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
| **Paralelismo con K workers** | Cada proceso tiene su propio caché: es seguro repartir la rejilla en K bloques contiguos (K = núcleos − 1), cada uno con su `worker_lote.py`. Speed-up ≈ K en rejillas grandes | Medio |
| **Orden serpenteante** | `itertools.product` salta del último valor de una fila al primero de la siguiente: la semilla h₁ queda lejos justo ahí. Recorrer en *boustrophedon* (invertir el sentido de la variable interna en cada fila) mantiene todos los vecinos adyacentes | Bajo |
| **Semilla del vecino, no del anterior** | En rejillas 2D+ el caso previo en la lista no siempre es el vecino más próximo. Guardar h₁ por índice de rejilla y sembrar con el vecino ya resuelto más cercano | Medio |
| **Riesgo de bloqueo por `stderr`** | El worker se lanza con `stderr=PIPE` pero `stderr` no se lee hasta que el proceso termina. Si escribe más que el búfer del pipe (≈ 4–64 KB; p. ej. muchos `RuntimeWarning` de numpy/scipy), el worker se bloquea y el caso termina en un **TIEMPO_AGOTADO falso**. Leer stderr en otro hilo o redirigirlo a un archivo | Bajo |
| **Excel incremental** | `guardar_en_excel` es O(filas totales) por guardado y crece sin fin. Usar `openpyxl` en modo *append* sobre la hoja, o un CSV/Parquet de trabajo que se exporte a Excel bajo demanda | Bajo |
| **Relanzar tras crash** | Tras un fin inesperado se marca ERROR solo el caso en curso y se relanza; si el crash es sistemático (p. ej. importación fallida) se relanza N veces. Abortar si el worker muere sin producir ninguna línea dos veces seguidas | Bajo |
| **`pendientes.remove(i)`** | O(n) por caso, O(n²) total; irrelevante hasta miles de casos, pero un `set` o un índice lo resuelve | Trivial |
| **Refinamiento de FRONTERA** | La frontera se localiza con la resolución de la rejilla. Una segunda pasada automática que biseque solo los tramos con cruce (sin interpolar, resolviendo casos reales) da la frontera a la precisión deseada con log₂ de las evaluaciones | Medio |
| **`worker_caso.py`** | Ya no lo usa la interfaz (`resolver_caso` no se llama desde `app.py`). Mantenerlo solo si se quiere el aislamiento total por caso; si no, eliminarlo evita duplicar lógica con `worker_lote.py` | Trivial |
