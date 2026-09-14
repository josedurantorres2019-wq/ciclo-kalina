# 03 — Inversión h→T y s→T

**Archivo:** [`kalina.py`](../kalina.py) — función [`estado_de`](../kalina.py#L161)

## Qué resuelve

El ciclo casi nunca conoce T directamente: conoce h (tras un balance de energía)
o s (tras una expansión isentrópica). `estado_de(P, w, h=…)` o
`estado_de(P, w, s=…)` encuentra la T que da ese valor.

## Física aplicada

A presión y composición fijas, h y s son **estrictamente crecientes** en T,
también dentro de la campana bifásica:

$$
\left(\frac{\partial h}{\partial T}\right)_{P,w} = c_p > 0,
\qquad
\left(\frac{\partial s}{\partial T}\right)_{P,w} = \frac{c_p}{T} > 0
$$

En la mezcla zeotrópica el cambio de fase ocurre con **deslizamiento de
temperatura** (entre T_burbuja y T_rocío), así que h(T) es continua pero con
quiebres de pendiente en esos dos puntos. La monotonía garantiza una raíz única;
los quiebres desaconsejan métodos que asumen derivada continua sin salvaguarda.

## Algoritmo

Busca la raíz de f(T) = estado(T,P,w)[h o s] − valor:

1. **Intervalo inicial:**
   - si ya se resolvió esa misma (P, w, propiedad): `[T₀ − 8, T₀ + 8]` K alrededor
     de la última T resuelta (memoria `_ultimaT`);
   - si no: `[max(235, T_sat,NH₃ − 60), T_sat,H₂O + 40]`.
2. **Ensanchamiento geométrico:** si f(lo) > 0 baja lo en pasos 16·1.6ᵏ K (piso
   233 K); si f(hi) < 0 sube hi igual (techo 900 K). Máx. 60 pasos en total.
3. **Brent** (`brentq`, `xtol = TOL_T = 1e-6 K`).
4. Guarda T en `_ultimaT` y devuelve el estado completo (**una evaluación extra**).

### Convergencia de Brent

Combina bisección, secante e interpolación cuadrática inversa. Garantiza no ser
peor que la bisección:

$$
n_{\max} \approx \log_2\frac{hi-lo}{x_{tol}} = \log_2\frac{16}{10^{-6}} \approx 24
$$

y en funciones suaves converge de forma superlineal (típicamente 6–10
evaluaciones con el intervalo tibio de 16 K).

### Piso de ruido

`TOL_T = 1e-6 K` con c_p ≈ 4 kJ/kg·K fija un piso de ≈ 4·10⁻⁶ kJ/kg para
cualquier h. Las tolerancias de los lazos superiores (`TOL_FRIO = 2e-5`,
`TOL_H1 = 1e-4`) están deliberadamente **por encima** de ese piso (correcciones
C3, C23bis del proyecto): pedir más precisión no converge, solo gasta
iteraciones.

## Costo

Cada evaluación de f es un `estado()` completo (doc 02). Con un flash frío y una
composición nueva, el código documenta ≈ 2–4 s por inversión; es el cuello de
botella de la rama Kalina (correcciones C27 del proyecto).

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Newton salvaguardado con c_p** | En fase única `prop` ya devuelve `cp`. Paso T ← T − (h − h*)/c_p, con respaldo a bisección si sale del intervalo (método `rtsafe`). En líquido/vapor converge en 2–3 evaluaciones | Medio |
| **Paso de secante desde la memoria** | Guardar en `_ultimaT` también el valor h(T₀): con dos puntos previos la primera estimación ya es secante, y el intervalo puede ser de ±1 K en lugar de ±8 K | Bajo |
| **Evitar la evaluación final** | `brentq` ya evaluó `estado(T)` en la raíz; envolver f para guardar el último dict evaluado y devolverlo en lugar de recalcularlo | Bajo |
| **Evaluaciones de acotamiento** | Siempre evalúa f(lo) y f(hi) aunque la memoria sea buena. Evaluar primero f(T₀) y expandir solo hacia el lado necesario ahorra 1 evaluación | Bajo |
| **Clave redondeada `round(w, 8)`** | En la sonda de incertidumbre y en barridos de composición cada w nuevo pierde la memoria. Usar la T de la w más cercana ya resuelta como semilla | Bajo |
| **Tabla h(T) por isobara** | Para (P, w) muy repetidos (estados 1, 10 del lazo frío), tabular h(T) una vez y usar interpolación como semilla de Newton | Medio |
