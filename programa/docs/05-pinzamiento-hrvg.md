# 05 — Pinzamiento en el HRVG con fuente de capacidad finita

**Archivo:** [`kalina.py`](../kalina.py) — funciones [`perfil_HRVG`](../kalina.py#L276), [`pinzamiento`](../kalina.py#L295), [`C_g_minimo_HRVG`](../kalina.py#L339)
**Se activa con:** `par['C_g']` (y opcionalmente `par['dT_pp_gas']`)

## Qué resuelve

Con la fuente como reservorio isotermo (hipótesis H10, `C_g` no declarado) basta
comprobar los extremos del HRVG. Con un gas real de capacidad calorífica finita,
las curvas pueden **cruzarse en el interior** del equipo aunque los extremos
estén bien. Este algoritmo detecta ese cruce y recorta el calor transferido hasta
que el equipo sea físicamente posible.

## Física aplicada

### Por qué ε-NTU no basta

En ε-NTU, Q_max = C_min (T_h,in − T_c,in). Si φ = C_g/C_wf ≥ 1, C_min es el
fluido de trabajo y Q_max no cambia con C_g: η sale idéntica para φ = 1, 2 o ∞.
Pero en una mezcla zeotrópica la curva h(T) del fluido tiene tres tramos de
pendiente muy distinta (líquido, evaporación con deslizamiento, vapor), mientras
que la del gas es una **recta**. Lo que limita es la **forma** del perfil, no la
capacidad total.

### Perfil en contracorriente

El gas entra caliente por el extremo del estado 2 y sale por el del estado 1.
Por balance de energía entre la sección h y el extremo caliente:

$$
T_g(h) = T_f - \frac{h_2 - h}{C_g},\qquad
T_{g,out} = T_f - \frac{Q_i}{C_g}
$$

**Segunda ley en todo punto** (con margen de diseño ΔT_pp):

$$
\Delta T(h) = T_g(h) - T_{wf}(h) \;\ge\; \Delta T_{pp}\qquad \forall\, h\in[h_1,h_2]
$$

### Cota de h₂ por nodo

Despejando h₂ para un nodo interior i fijo:

$$
T_f - \frac{h_2 - h_i}{C_g} - T_i \ge \Delta T_{pp}
\;\Longleftrightarrow\;
h_2 \le b_i = h_i + C_g\,(T_f - T_i - \Delta T_{pp})
$$

### Capacidad mínima de la fuente

Con el ciclo ya resuelto, el menor C_g que no cruza:

$$
C_{g,\min} = \max_i \frac{h_2 - h_i}{T_f - T_i - \Delta T_{pp}}
$$

## Algoritmo

1. **Tabular el perfil una sola vez** en N = 40 nodos equiespaciados en
   **temperatura** entre T₁ y T₂: h(T) es evaluación directa, T(h) exigiría una
   inversión por nodo. Los extremos se fuerzan a los valores reales (T₁,h₁),
   (T₂,h₂) para no perder el mínimo del extremo caliente por un artefacto de
   rejilla.
2. Calcular ΔT en todos los nodos; si min ΔT ≥ ΔT_pp, no hay recorte.
3. Si lo hay: calcular bᵢ y su **mínimo acumulado** desde el extremo frío,
   mᵢ = min(b₀…bᵢ). Un h₂ es factible mientras h₂ ≤ mᵢ para todos los nodos que
   quedan dentro del equipo (hᵢ ≤ h₂).
4. El primer nodo k con hₖ > mₖ marca dónde deja de ser factible; el h₂ recortado
   se obtiene por **interpolación lineal** entre k−1 y k sobre la misma tabla.
5. `resolver` sustituye el estado 2 por `estado(T2_recortado)` y sigue.

Complejidad: O(N) evaluaciones de `estado` + O(N) aritmética vectorizada con
numpy. No hay lazo de Brent interno.

## Costo

40 evaluaciones de `estado()` **por cada evaluación de `cierre_hrvg`**, es decir,
por cada evaluación de `tramo(h1)` del lazo exterior (doc 04). Suele dominar el
costo cuando `C_g` está declarado.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Rejilla adaptativa** | El mínimo de ΔT aparece en los quiebres de h(T) (T_burbuja, T_rocío). Calcular esas dos T una vez y concentrar nodos ahí; en tramos monofásicos h(T) es casi lineal y bastan 3–4 nodos. N efectivo ≈ 12–15 | Medio |
| **Reutilizar el perfil entre evaluaciones de h₁** | Para la misma (P_alta, w_b) h(T) no depende de h₁: tabular una vez sobre [T_min, T_f] y extraer el tramo [T₁, T₂] por interpolación. Convierte 40·n_tramo evaluaciones en ≈ 40 por caso | Bajo–medio |
| **Nodos ya calculados en `_flash`** | Con T equiespaciado y P fija, los flashes quedan en la tabla del doc 02; pero `_mono` sigue recalculando densidades. Cachear h(T,P,w) completo por nodo redondeado | Bajo |
| **Refinamiento local tras el recorte** | La interpolación lineal entre nodos introduce un error O(ΔT²·h''); un paso de secante con 1–2 evaluaciones reales alrededor del h₂ recortado lo elimina sin subir N | Bajo |
| **`C_g_minimo_HRVG` sin perfil extra** | Si el limitador ya tabuló el perfil en la última evaluación, reutilizarlo en lugar de volver a llamar `perfil_HRVG` al final de `resolver` | Bajo |
