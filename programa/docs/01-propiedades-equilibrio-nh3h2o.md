# 01 — Propiedades y equilibrio líquido-vapor NH₃-H₂O

**Archivo:** [`nh3h2o.py`](../nh3h2o.py) · **Capa:** 0 (todo lo demás se apoya aquí)

## Qué resuelve

Dado (T, ρ, x) calcula todas las propiedades de la mezcla; dado (T, P, x)
encuentra la densidad; y dado (T, x) encuentra la presión de burbuja o de rocío
y la composición de la fase conjugada. `x` es fracción **molar** de NH₃.

## Física aplicada

### Ecuación fundamental (IAPWS G4-01)

La mezcla se describe con la energía de Helmholtz adimensional, suma de una
parte ideal y una residual:

$$
\frac{a(T,\rho,x)}{RT} = \Phi(\tau,\delta,x) = \Phi^{\circ}(\tau,\delta,x) + \Phi^{r}(\tau,\delta,x),
\qquad \tau=\frac{T_n(x)}{T},\quad \delta=\frac{\rho}{\rho_n(x)}
$$

con la parte residual construida a partir de las de los puros más un término
de desviación:

$$
\Phi^{r} = (1-x)\,\Phi^{r}_{\mathrm{H_2O}} + x\,\Phi^{r}_{\mathrm{NH_3}} + \Delta\Phi^{r},
\qquad
\Delta\Phi^{r} = f(x)\left[A(\tau,\delta) + x\,B(\tau,\delta) + x^{2}C(\tau,\delta)\right]
$$

Todas las propiedades salen de derivadas de Φ. La presión, por ejemplo
([`P_of`](../nh3h2o.py#L83)):

$$
P = \rho R T \left(1 + \delta\,\frac{\partial \Phi^{r}}{\partial \delta}\right) = Z\,\rho R T
$$

### Equilibrio de fases

Dos fases coexisten a (T, P) cuando la fugacidad de **cada componente** es igual
en ambas:

$$
x_i\,\varphi_i^{L}(T,P,x) = y_i\,\varphi_i^{V}(T,P,y),\qquad i\in\{\mathrm{NH_3},\mathrm{H_2O}\}
$$

Los coeficientes de fugacidad dependen de la derivada composicional
∂Φʳ/∂x. Por eso un error en esa derivada contamina **todo** el equilibrio
aunque las propiedades monofásicas sean correctas.

### Correcciones al paquete `iapws` 1.5.5

[`_aplicar_correcciones`](../nh3h2o.py#L48) reescribe dos métodos del paquete al
importarse (sin tocar la instalación):

- **P1** — errata en el exponente α de la Ec. (4): `1.12455` → `1.125455`.
- **P2** — faltaba la regla del producto en la derivada de la Ec. (8):

$$
\frac{\partial \Delta\Phi^{r}}{\partial x} = f'(x)\left[A + xB + x^{2}C\right] + f(x)\left[B + 2xC\right]
$$

  El paquete solo tenía el primer sumando. Esto explicaba desviaciones de 0.5–4 %
  en presiones de burbuja y rocío.

Los parches se aplican por reemplazo exacto de texto y **abortan** si el texto no
aparece exactamente una vez: por eso `iapws` debe ser 1.5.5.

## Algoritmos

### A. Raíz de densidad — [`rho_TPx`](../nh3h2o.py#L89)

Resuelve P(ρ) − P = 0 con `brentq` después de acotar la raíz:

- **Vapor:** parte de la densidad de gas ideal ρ₀ = 0.98·P/(RT) y multiplica
  por 1.06 hasta que cambia el signo (máx. 200 pasos).
- **Líquido:** baja desde 70 mol/dm³ en pasos de 2 hasta un valor finito con
  g ≥ 0, y luego multiplica por 0.97 hasta g < 0 (máx. 300 pasos).

Cada paso de acotamiento evalúa Φʳ completo.

### B. Burbuja y rocío — [`_flashP`](../nh3h2o.py#L125) + [`_polish`](../nh3h2o.py#L162)

1. **Inicialización de Raoult** con las presiones de saturación de los puros:

$$
P^{(0)}_{\text{bub}} = z\,p_A^{sat} + (1-z)\,p_W^{sat},\qquad
P^{(0)}_{\text{dew}} = \left(\frac{z}{p_A^{sat}}+\frac{1-z}{p_W^{sat}}\right)^{-1}
$$

2. **Sustitución sucesiva** con valores K = φᴸ/φⱽ (para burbuja):

$$
K_i=\frac{\varphi_i^{L}}{\varphi_i^{V}},\qquad
S=\sum_i z_i K_i,\qquad
P^{(k+1)} = P^{(k)}\,S,\qquad
y_i^{(k+1)}=\frac{z_i K_i}{S}
$$

   con relajación 0.6/0.4 y criterio |ΔP|/P + |Δy| < 10⁻¹⁰.

3. **Pulido de Newton** (`fsolve`) sobre el sistema exacto en (ρ_L, ρ_V, y):
   igualdad de presiones y de fugacidades. Necesario cerca del locus crítico,
   donde la sustitución sucesiva se estanca.

### C. Temperatura de saturación — [`_satT`](../nh3h2o.py#L201)

`brentq` sobre T de `bubbleP(T,x) − P`, acotado por las temperaturas de
saturación de los puros (la mezcla no tiene azeótropo, así que
T_sat,NH₃ < T_bub, T_dew < T_sat,H₂O).

**Anidamiento:** cada evaluación de `bubbleT` = ~15 iteraciones de Brent × (hasta
300 pasos de sustitución + un `fsolve`) × 2 raíces de densidad × acotamiento. Es
la operación más cara del módulo.

## Costo y dónde optimizar

| Oportunidad | Por qué | Esfuerzo |
|---|---|---|
| **Newton en `rho_TPx`** con ∂P/∂ρ = RT(1 + 2δΦʳ_δ + δ²Φʳ_δδ) | `_phir` ya calcula `fird`, `firdd`; convergencia cuadrática frente al acotamiento geométrico (decenas de evaluaciones de Φʳ por raíz) | Medio |
| **Semilla de densidad** memorizada por (T, P, x, fase) redondeados | En los lazos del ciclo se piden densidades casi idénticas miles de veces | Bajo |
| **Aceleración de la sustitución sucesiva** (GDEM o Newton completo desde el inicio) | La relajación fija 0.6/0.4 converge linealmente; cerca de la zona crítica, muy lento | Medio |
| **Jacobiano analítico en `_polish`** | `fsolve` estima el Jacobiano 3×3 con diferencias finitas: 3 evaluaciones extra de `prop` (cada una con 2 raíces de densidad) por iteración | Alto |
| **Límites de caché** en `_psat_cache` | Crece sin límite en barridos largos; `functools.lru_cache(maxsize=…)` | Bajo |
| **Vectorizar Φʳ** sobre arreglos de ρ | Los sumatorios de `iapws` son bucles de Python escalares | Alto (reescribir la evaluación) |

> Cualquier optimización debe volver a pasar `validacion_g4_01.py` sin cambiar
> las desviaciones reportadas.
