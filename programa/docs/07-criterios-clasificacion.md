# 07 — Criterios de admisibilidad y clasificación

**Archivos:** [`kalina.py`](../kalina.py) — [`criterios`](../kalina.py#L756), [`clasificar`](../kalina.py#L896), [`supuestos`](../kalina.py#L944); [`sugerencias.py`](../sugerencias.py)

## Qué resuelve

Convierte una solución numérica en un **diagnóstico**: qué leyes físicas,
restricciones técnicas o condiciones topológicas cumple o incumple el punto, qué
etiqueta le corresponde y qué incertidumbre dejan las variables no declaradas.
Es el entregable principal del programa.

## Física aplicada (criterios)

| Id | Severidad | Condición | Fundamento |
|---|---|---|---|
| N2 | numérico | \|Q_i − Q_out − W_net\| ≤ 10⁻⁶ Q_i | 1.ª ley global |
| N3 | numérico | \|w_b − m_r w_r − m_p w_p\| ≤ 10⁻⁶ w_b | Conservación de NH₃ |
| S1 | crítico | s₄ ≥ s₃ | 2.ª ley en la turbina adiabática |
| S2 | crítico | s₁₀ ≥ s₉ | 2.ª ley en la bomba |
| S3 | crítico | s₇ ≥ s₆ | 2.ª ley en la válvula (no aplica si degenerado) |
| S4 | crítico | T₅ > T₁ y T₆ > T₁₀ | Calor fluye de caliente a frío en el regenerador |
| S5 | crítico | T₂ ≤ T_f | Extremo caliente del HRVG |
| S6 | crítico | T₉ ≥ T_agua,in | Extremo frío del condensador |
| S8 | crítico | η < η_Carnot | 2.ª ley global |
| S9 | crítico | min ΔT_HRVG ≥ ΔT_pp,gas | Sin cruce interno (doc 05) |
| S10 | crítico | min ΔT_cond ≥ ΔT_pp,cond | Sin cruce interno (doc 06) |
| O1 | tecnológico | q₄ ≥ 0.90 | Erosión de álabes |
| O2 | crítico | estado 9 líquido | Cavitación de la bomba |
| O3 | crítico | W_net > 0 | Ciclo motor |
| O5 | topológico | 10⁻³ ≤ q₂ ≤ 1 − 10⁻³ | El separador separa |
| O6 | crítico | T_gas,out ≥ T_min,gas | Rocío ácido |
| O7 | crítico | T_agua,out ≤ T_max | Vertido |
| C1 | crítico | w_p < w_b < w_r | Orden de composiciones |
| C1b | advertencia | w_r − w_p > 0.01 | Separación mayor que la banda de G4-01 |
| C3 | crítico | 0 ≤ m_r, m_p ≤ 1 | Caudales físicos |
| F1, PP | advertencia | T₁ < T_sat,L; pinch en extremo caliente | Coherencia del cierre del regenerador |
| CD, CDc | advertencia | cierre declarado alcanzado | El limitador recortó h₂ o h₉ |

Fórmulas usadas por los criterios:

$$
\eta_{Carnot} = 1 - \frac{T_{amb}}{T_f},\qquad
\Delta S_{turbina} = s_4 - s_3 \ge 0,\qquad
Q_i - Q_{out} - W_{net} = 0
$$

## Algoritmo de clasificación

Jerarquía por severidad del criterio incumplido más grave:

```
numerico    → NO_CONVERGE
critico     → INVIABLE
topologico  → DEGENERADO
tecnologico / advertencia → VALIDO con advertencia
ninguno     → VALIDO
```

La lista de violaciones se conserva completa: un punto puede ser a la vez
DEGENERADO e INVIABLE, y ambos hallazgos se reportan.

Los criterios que dependen de variables opcionales **solo se agregan** si la
variable está declarada; la ausencia se documenta en `supuestos()`, que para cada
variable opcional devuelve el valor asumido y la incertidumbre medida que deja.

`sugerencias.py` es un diccionario criterio → texto accionable, aplicado sobre la
lista de violados.

## Costo

Despreciable frente a `resolver()`, salvo dos llamadas al motor:

- `T_sat(P_alta, w_b, 'L')` en F1 (una `bubbleT`, cacheada en `_env`).
- `diagnostico_condensador()` en `supuestos()` si `calcular_condensador=True`
  (un perfil completo del condensador, doc 06).

## Dónde optimizar / mejorar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Criterios como datos** | Cada criterio está escrito como código imperativo. Una tabla `(id, severidad, predicado, formato, aplica_si)` permitiría que la interfaz liste todos, que `sugerencias.py` compruebe cobertura y que se prueben uno por uno | Medio |
| **Tolerancias con nombre** | `1e-6`, `1e-9`, `0.90`, `0.01` aparecen literales; subirlas a constantes junto a `TOL_T` hace explícito su vínculo con el piso de ruido | Bajo |
| **N2 relativo al piso de ruido** | N2 exige 10⁻⁶·Q_i; con Q_i ≈ 1000 kJ/kg eso son 10⁻³ kJ/kg, holgado frente a TOL_H1 = 10⁻⁴, pero conviene derivarlo de TOL_H1 para que no quede desfasado si se cambian las tolerancias | Bajo |
| **F1 sin `bubbleT`** | La T de burbuja a P_alta y w_b puede obtenerse del perfil del HRVG (quiebre de h(T)) cuando está tabulado, evitando un `bubbleT` completo | Bajo |
| **Cobertura de `sugerencias.py`** | Faltan S1–S8, N2, N3, O7, F1, PP, CD, CDc; decisión consciente del proyecto (ampliar con evidencia), pero una prueba que liste los no cubiertos ayuda a seguirlo | Bajo |
