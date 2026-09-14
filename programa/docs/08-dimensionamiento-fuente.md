# 08 — Dimensionamiento a la fuente

**Archivo:** [`kalina.py`](../kalina.py) — función [`dimensionar`](../kalina.py#L1039)

## Qué resuelve

Dada una fuente real (capacidad calorífica absoluta C_g,abs en kW/K y una
temperatura mínima de salida del gas por rocío ácido), encuentra el **mayor
caudal de solución básica m_b** que la planta puede tener.

## Física aplicada

### m_b no es un grado de libertad termodinámico

Todas las ecuaciones del ciclo son específicas (por kg/s de m_b). m_b solo
multiplica Q_i, W_t, W_p, W_net. Lo que cambia la física es la **razón**:

$$
C_g = \frac{C_{g,abs}}{\dot m_b}\quad\text{[kW/K por kg/s de } m_b\text{]}
$$

Aumentar m_b equivale a una fuente relativamente más pobre.

### Restricción 1 — rocío ácido

$$
T_{g,out}(\dot m_b) = T_f - \frac{Q_i}{C_{g,abs}/\dot m_b} \;\ge\; T_{min,gas}
$$

T_g,out decrece con m_b: más caudal extrae más calor del mismo gas.

### Restricción 2 — cierre declarado sin recorte

$$
\Delta T_{min,HRVG}(\dot m_b) \ge \Delta T_{pp,gas}
$$

El mayor m_b para el que el limitador del doc 05 aún no recorta h₂, es decir, el
cierre declarado (ΔT_app o ε_HRVG) sigue describiendo el equipo.

### Resultado

$$
\dot m_b^{*} = \min\big(\dot m_{b,rocio},\;\dot m_{b,pinch}\big)
$$

y `cual_ato` indica qué restricción manda.

## Algoritmo

Para cada restricción g(m_b) = 0 ([`raiz`](../kalina.py#L1080)):

1. Evalúa g en `mb_lo = 0.05` y `mb_hi = max(C_g,abs, 1)`.
2. Si no hay cambio de signo, multiplica `hi` por 1.5 hasta 30 veces.
3. Si sigue sin cambio de signo, devuelve el extremo (la restricción no ata en el
   rango); si lo hay, `brentq(xtol = 1e-4 kg/s)`.
4. Una evaluación final de `resolver()` en m_b*.

Cada evaluación de g es un `resolver()` completo con `C_g` declarado.

## Costo

(2 + expansiones + iteraciones de Brent) × 2 restricciones + 1, **cada una un
`resolver()` completo** con el limitador activo. Del orden de 20–40 resoluciones
del ciclo.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Memoizar `evalua(mb)`** | `g_rocio` y `g_pinch` llaman a `resolver()` en los **mismos** m_b iniciales (0.05, mb_hi y las expansiones): se resuelven dos veces. Un dict m_b → (est, ind) elimina ≈ 1/3 de las resoluciones, y la evaluación final en m_b* casi siempre ya está calculada | Bajo |
| **Encadenar `h1_semilla`** | Cada `resolver()` parte del rango completo de h₁. Pasar el h₁ de la evaluación anterior (como hace `worker_lote.py`) acelera ~15–20× según `incertidumbre.py` | Bajo |
| **Una sola búsqueda** | Buscar la raíz de min(g_rocio, g_pinch) (función continua) resuelve ambas a la vez, y luego se identifica cuál es activa en m_b* | Bajo |
| **Tolerancia relativa** | `xtol = 1e-4 kg/s` absoluto es excesivo para plantas de decenas de kg/s y escaso para m_b ≈ 0.05; usar `rtol` | Bajo |
| **Mejor semilla analítica para rocío** | Si Q_i varía poco con m_b (fuera del régimen de recorte), m_b,rocío ≈ C_g,abs (T_f − T_min)/Q_i da un intervalo de ±10 % en lugar de expandir desde 0.05 | Bajo |
| **Interfaz** | `dimensionar` no está expuesta en `app.py`; hoy solo se usa desde scripts | Medio |
