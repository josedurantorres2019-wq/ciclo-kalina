# 04 — Resolución del ciclo KCS-11

**Archivo:** [`kalina.py`](../kalina.py) — función [`resolver`](../kalina.py#L457)

## Qué resuelve

Encuentra los 10 estados termodinámicos del ciclo y sus indicadores (Q_i, W_net,
η, m_r, …) para un conjunto de parámetros `par`, con caudal de solución básica
normalizado a m_b = 1.

## El ciclo

```
          fuente (T_f)
              │
 1 ──► [HRVG] ──► 2 ──► [Separador] ──► 3 (rica, vapor) ──► [Turbina] ──► 4 ─┐
 ▲                          │                                                 │
 │                          ▼ 5 (pobre, líquido)                              │
 │                    [Regenerador] ──► 6 ──► [Válvula] ──► 7 ──► [Absorbedor]◄┘
 │                          ▲                                        │
 │                          │                                        ▼ 8
 └──────────── 10 ◄── [Bomba] ◄── 9 ◄──────────── [Condensador] ◄────┘
                                                   sumidero (T_amb)
```

La solución básica fría (10) se precalienta en el regenerador con la pobre
caliente (5 → 6) antes de entrar al HRVG como estado 1.

## Física aplicada (balances por equipo)

Todos los equipos en régimen permanente, sin pérdidas de calor ni de presión.

**Separador** (estado 2 bifásico; conservación de masa total y de NH₃):

$$
w_r = w_V(T_2,P_{alta}),\quad w_p = w_L(T_2,P_{alta}),\qquad
m_r = \frac{w_b - w_p}{w_r - w_p},\quad m_p = 1 - m_r
$$

Si el estado 2 **no** es bifásico, el separador no separa: m_p = 0, w_r = w_p = w_b
(**rama degenerada**, mismo camino de código).

**Turbina** (eficiencia isentrópica):

$$
s_{4s} = s_3,\qquad h_4 = h_3 - \eta_t\,(h_3 - h_{4s})
$$

**Bomba:**

$$
s_{10s} = s_9,\qquad h_{10} = h_9 + \frac{h_{10s} - h_9}{\eta_p}
$$

**Válvula** (isentálpica): h₇ = h₆ a P_baja.

**Absorbedor** (mezcla adiabática): h₈ = m_r h₄ + m_p h₇.

**Regenerador** (lado frío = básica, lado caliente = pobre):

$$
h_1 = h_{10} + m_p\,(h_5 - h_6)
$$

**Cierres** (grados de libertad de cada intercambiador; el usuario elige uno):

| Equipo | Método A (diferencia de T) | Método B (efectividad) |
|---|---|---|
| HRVG | T₂ = T_f − ΔT_app | h₂ = h₁ + ε_H (h₂,ideal − h₁), h₂,ideal = h(T_f) |
| Regenerador | T₆ = T₁₀ + ΔT_pp | h₆ = h₅ − ε_R (h₅ − h₆,ideal), h₆,ideal = h(T₁₀) |
| Condensador | T₉ = T_agua,in + ΔT_app | h₉ = h₈ − ε_C (h₈ − h₉,ideal), h₉,ideal = h(T_agua,in) |

**Indicadores:**

$$
Q_i = h_2 - h_1,\quad Q_{out} = h_8 - h_9,\quad
W_t = m_r(h_3 - h_4),\quad W_p = h_{10} - h_9,\quad W_{net} = W_t - W_p
$$

$$
\eta = \frac{W_{net}}{Q_i},\qquad \eta_{Carnot} = 1 - \frac{T_{amb}}{T_f},\qquad
\eta_{II} = \frac{\eta}{\eta_{Carnot}},\qquad RBW = \frac{W_p}{W_t}
$$

## Algoritmo: dos lazos anidados

El ciclo tiene dos corrientes de "corte" (variables desconocidas que cierran un
recorrido cerrado):

### Lazo exterior — punto fijo en h₁ resuelto con Brent

La función [`tramo(h1)`](../kalina.py#L605) recorre el ciclo desde h₁ y devuelve
el h₁ que resulta del regenerador, G(h₁). Se resuelve

$$
F(h_1) = G(h_1) - h_1 = 0,\qquad h_1 \in \big[h(T_{agua,in}),\; h_{2,ideal}\big]
$$

con `brentq(xtol = TOL_H1 = 1e-4 kJ/kg)`.

- La ausencia de cambio de signo en el intervalo **es** la definición operativa
  de `NO_CONVERGE`.
- **Método A:** T₂ es dato, G no depende de h₁. El resolvedor lo detecta porque
  G(lo) = G(hi) y devuelve sin iterar.
- **Arranque tibio:** si `par['h1_semilla']` existe, el intervalo inicial es
  semilla ± 2 kJ/kg, ensanchado ×3 hasta tener cambio de signo (máx. 15 veces).

### Lazo interior frío — sustitución sucesiva en h₁₀

Dentro de cada `tramo`, [`lazo_frio`](../kalina.py#L555) resuelve
h₁₀ → T₁₀ → h₆ → h₇ → h₈ → h₉ → h₁₀. Cuando el condensador cierra por efectividad,
h₉ depende de h₈, que depende de h₆, que depende de T₁₀: aparece un bucle. Su
ganancia (según el código):

$$
g \approx (1-\varepsilon_C)\,\varepsilon_R\,m_p \approx 0.14
\quad\Rightarrow\quad
|e_k| \approx g^{k}|e_0|,\qquad
k \approx \frac{\ln(\text{tol}/e_0)}{\ln g}
$$

Con e₀ ≈ 5 kJ/kg y tol = 2·10⁻⁵ salen ≈ 6–10 pasadas. Con el condensador por
approach, g = 0 y termina en una pasada.

Salidas del lazo: residuo < `TOL_FRIO` (2e-5), o estancado 6 pasadas con residuo
< `TOL_FRIO_PISO` (1e-4). Después recalcula los estados con el h₁₀ final.

## Costo

Por evaluación de `tramo`: ≈ 3 inversiones (T₁, e4s, e4) + pasadas del lazo frío
× 4 inversiones (e10, e6 o e7, e9, e10s) + los perfiles de pinzamiento si están
activos (40 evaluaciones de `estado` por pasada del HRVG, ver doc 05).

Brent sobre h₁ típicamente requiere 8–12 evaluaciones de `tramo`, más las 2 de
acotamiento, más 1 final. Todo el costo es multiplicativo.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Aitken Δ² / Wegstein en el lazo frío** | Con ganancia ≈ 0.14 constante, la extrapolación de Aitken h* ≈ h_k − (Δh_k)²/(Δ²h_k) reduce ~8 pasadas a 3–4; cada pasada son 4 inversiones | Bajo |
| **Semilla de h₁₀ entre evaluaciones de `tramo`** | Hoy parte siempre de h₉,ideal + 5. Guardar el h₁₀ convergido de la llamada anterior (Brent evalúa h₁ cercanos) deja e₀ ≈ 10⁻² y 1–2 pasadas | Bajo |
| **Resolver el lazo frío con Brent/secante** en vez de sustitución | Útil si ε_C baja y la ganancia se acerca a 1 (convergencia lenta) | Medio |
| **Recalcular al salir** | Tras converger, el lazo re-evalúa e10, e6, e7, e9 (≈ 4 inversiones) aunque el cambio sea < tol; se puede devolver el último juego calculado | Bajo |
| **`tramo(h1)` final repetido** | Tras `brentq` se llama de nuevo `tramo(h1)`; memoizar la última evaluación del Brent | Bajo |
| **Detección del método A antes de iterar** | Si `'dT_app' in cierre_HRVG` y no hay C_g, G es constante por construcción: se ahorra una evaluación de `tramo` | Bajo |
| **Newton/secante sobre h₁ con semilla** | Con semilla del punto vecino, la secante converge en 2–3 evaluaciones frente a las ≥ 4 de acotar + Brent | Medio |

> `prueba5.py` es la prueba de regresión de este algoritmo: resuelve el mismo
> punto por A (sin iterar) y por B (con ambos lazos) y exige coincidencia de
> 10⁻³ K / 10⁻³ kJ/kg. Cualquier cambio en los lazos debe pasarla.
