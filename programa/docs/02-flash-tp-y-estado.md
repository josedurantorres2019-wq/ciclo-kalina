# 02 — Flash a (T, P) y estado de la mezcla

**Archivo:** [`kalina.py`](../kalina.py) — funciones [`flash_TP`](../kalina.py#L34), [`estado`](../kalina.py#L95), [`Tsat_puros`](../kalina.py#L26)

## Qué resuelve

`estado(T, P, w)` devuelve h, s, título q y fase de la mezcla con composición
global **másica** w. Es la función más llamada del programa: cada inversión,
cada perfil y cada paso de cada lazo pasa por aquí.

## Física aplicada

### Conversión másica ↔ molar

$$
x = \frac{w/M_{A}}{w/M_{A} + (1-w)/M_{W}},\qquad
w = \frac{x\,M_{A}}{x\,M_{A} + (1-x)\,M_{W}}
$$

### Diagrama T-x a presión fija

A una presión P la mezcla zeotrópica NH₃-H₂O tiene, para cada T entre las
saturaciones de los puros, una composición de líquido saturado x_L(T,P) y otra de
vapor saturado x_V(T,P). La composición global x decide la fase:

$$
x \le x_L \Rightarrow \text{líquido},\qquad
x \ge x_V \Rightarrow \text{vapor},\qquad
x_L < x < x_V \Rightarrow \text{bifásico}
$$

### Regla de la palanca (en base másica)

$$
q = \frac{w - w_L}{w_V - w_L},\qquad
h = q\,h_V + (1-q)\,h_L,\qquad
s = q\,s_V + (1-q)\,s_L
$$

El título es másico porque h y s están en kJ/kg.

## Algoritmo

### `flash_TP(T, P)`

Resuelve el sistema 2×2 de igualdad de fugacidades en (x_L, x_V) con T y P fijos:

$$
F(x_L,x_V)=
\begin{bmatrix}
\ln(x_L\varphi_A^{L}) - \ln(x_V\varphi_A^{V})\\
\ln\big((1-x_L)\varphi_W^{L}\big) - \ln\big((1-x_V)\varphi_W^{V}\big)
\end{bmatrix} = 0
$$

1. **Banda de guarda:** si T está a menos de 0.3 K de la saturación de un puro,
   devuelve `None`. Ahí la ventana bifásica se cierra (x_L → 1 o x_V → 0) y el
   Jacobiano se vuelve singular.
2. **Tabla por presión:** `_flash[P]` es una lista ordenada por T de soluciones
   ya obtenidas. Búsqueda binaria (`bisect`); si T ya está, se devuelve.
3. **Arranques**, en orden:
   - *tibio interpolado* (2b): si la tabla tiene soluciones a **ambos** lados de
     T, interpola x_L, x_V linealmente entre los dos vecinos — la campana
     zeotrópica es casi lineal en T a P fija, la semilla queda mucho más cerca
     y `fsolve` converge en ~1-2 iteraciones en vez de ~3;
   - *tibio*: solo un vecino disponible → usa el más cercano;
   - *respaldo*: Raoult, x_L = (P − p_W)/(p_A − p_W), x_V = x_L·p_A/P.
4. `fsolve(xtol=1e-12)`; se acepta solo si converge y 0 < x_L < x_V < 1.
5. Inserta la solución en la tabla (`bisect.insort`).

El arranque tibio baja las iteraciones de `fsolve` de ~9 a ~3; la interpolación
(2b, implementada 2026-09-16) la baja a ~1-2 y reduce los reintentos con Raoult.
En B_base (barrido con tabla grande): FRIO 46.8 → 37.7 s (-19 %), vecino
34.3 → 28.1 s (-18 %). Resultado en el PLAN_OPTIMIZACION.

### `estado(T, P, w)`

- Si el flash devuelve solución: clasifica con x_L, x_V y aplica la palanca.
- Si no (banda de guarda o no convergencia): decide la fase con `bubbleP`/`dewP`
  de `nh3h2o`, que están parametrizadas por x y sí son estables junto a los
  puros; si P queda entre las dos, lanza `RuntimeError`.

## Costo

Cada evaluación de F = 2 raíces de densidad + 2 llamadas a `prop` (ver doc 01).
`fsolve` sin Jacobiano añade 2 evaluaciones de F por iteración para estimarlo
por diferencias finitas. Un flash nuevo cuesta ≈ (3 iteraciones + Jacobiano) ×
4 densidades.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Jacobiano analítico** | ∂ln φᵢ/∂x se obtiene de las segundas derivadas de Φʳ. Pasar `fprime` a `fsolve` o escribir Newton 2×2 explícito elimina las evaluaciones por diferencias finitas (≈ 2 de cada 3) | Alto |
| **Reconstrucción de `Ts` en cada llamada** | [`Ts = [r[0] for r in tabla]`](../kalina.py#L52) es O(n) por llamada y la tabla crece con el barrido. Mantener una lista paralela de T ordenada o usar `sortedcontainers` | Bajo |
| **Interpolación en lugar de solo semilla** | Si hay soluciones a ambos lados de T, interpolar x_L, x_V linealmente da una semilla mucho más cercana (a menudo converge en 1 iteración) | Bajo | **HECHO (2b)**: B_base -19 % FRIO, -18 % vecino |
| **Tablas sin límite** | `_flash`, `_puros`, `_env`, `_ultimaT` crecen sin tope. En un proceso persistente (worker de lote) conviene un LRU por clave de presión | Bajo |
| **Densidades de la iteración anterior** | Dentro de `fsolve` las composiciones cambian poco: sembrar `rho_TPx` con la última densidad por fase | Medio |
| **Evaluación duplicada en `estado`** | En la rama bifásica `_mono(T,P,xL,'l')` y `_mono(T,P,xV,'v')` recalculan densidades que el propio flash ya obtuvo en su última evaluación de F; devolverlas desde el flash ahorra 2 raíces + 2 `prop` por llamada | Medio |

> **Cuidado con los cachés:** `resolver()` evita deliberadamente llamadas extra
> (p. ej. no calcula `C_cold_min`) porque llenar las tablas de arranque tibio
> mueve los resultados dentro del ruido numérico. Cualquier cambio de caché debe
> comprobarse con `prueba5.py`.
