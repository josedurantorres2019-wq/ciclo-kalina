# 06 — Condensador con agua de capacidad finita

**Archivo:** [`kalina.py`](../kalina.py) — funciones [`perfil_condensador`](../kalina.py#L367), [`pinzamiento_condensador`](../kalina.py#L402), [`C_cold_minimo`](../kalina.py#L441), [`diagnostico_condensador`](../kalina.py#L929)
**Se activa con:** `par['C_cold']` (y opcionalmente `T_agua_in`, `dT_pp_cond`)

## Qué resuelve

Es el espejo del doc 05 en el lado frío. Con sumidero infinito solo importa el
extremo frío (T₉). Con un caudal de agua finito, el agua se calienta a lo largo
del equipo y la mezcla, que condensa con deslizamiento, puede acercarse al agua
**en el interior**. Si el margen no se cumple, se retira menos calor (h₉ sube).

## Física aplicada

Contracorriente: el agua entra fría por el extremo del estado 9 y sale por el
del estado 8.

$$
T_w(h) = T_{w,in} + \frac{h - h_9}{C_{cold}},\qquad
T_{w,out} = T_{w,in} + \frac{Q_{out}}{C_{cold}},\qquad
\dot m_{agua} = \frac{C_{cold}}{c_{p,agua}}
$$

Segunda ley con margen:

$$
T_{mix}(h) - T_w(h) \ge \Delta T_{pp,cond}\qquad \forall\,h\in[h_9,h_8]
$$

Cota de h₉ por nodo (ahora es un **mínimo** para h₉):

$$
T_i - T_{w,in} - \frac{h_i - h_9}{C_{cold}} \ge \Delta T_{pp}
\;\Longleftrightarrow\;
h_9 \ge b_i = h_i - C_{cold}\,(T_i - T_{w,in} - \Delta T_{pp})
$$

Capacidad mínima de agua:

$$
C_{cold,\min} = \max_i \frac{h_i - h_9}{T_i - T_{w,in} - \Delta T_{pp}}
\quad(\infty \text{ si el propio extremo frío incumple})
$$

## Algoritmo

### Perfil (corrección C31)

La versión anterior invertía T₈ = T(h₈) en cada pasada del lazo frío, justo en la
zona donde el flash es frágil (cerca del rocío a P_baja), y el ciclo no
convergía. Ahora:

1. Sube en temperatura desde T₉ con paso fijo ΔT = 2 K (tramo máx. 300 K),
   evaluando h(T) directamente.
2. Si `estado()` falla en un nodo, lo **salta** y lo cuenta
   (`cond_nodos_saltados`, reportado en S10).
3. Descarta nodos no crecientes en h.
4. Al primer nodo con h ≥ h₈, interpola linealmente T₈ y termina: el último nodo
   es exactamente (T₈, h₈).

### Recorte

1. ΔTᵢ = Tᵢ − T_w(hᵢ). Si min ≥ ΔT_pp, no hay recorte.
2. bᵢ como arriba y su **máximo acumulado desde el extremo caliente**:
   Mᵢ = max(bᵢ…b_N).
3. gᵢ = hᵢ − Mᵢ es creciente; se busca su primer cruce por cero e interpola h₉.
4. Si ningún nodo cumple, `sin_solucion = True` (ningún h₉ satisface el margen
   con ese C_cold) y S10 lo reporta como crítico.

## Costo

El número de nodos depende del salto T₈ − T₉: ≈ (T₈ − T₉)/2 evaluaciones de
`estado()` (típicamente 20–60) **por cada pasada del lazo frío** × cada
evaluación de `tramo`. Es el término más caro cuando `C_cold` está declarado.

`C_cold_min` no se calcula dentro de `resolver()` (llenaría las memorias de
arranque tibio y movería resultados en el ruido); solo bajo demanda con la
casilla de la interfaz.

## Dónde optimizar

| Oportunidad | Detalle | Esfuerzo |
|---|---|---|
| **Perfil una vez por (P_baja, w_b)** | h(T) a P_baja y w_b no depende de las iteraciones: tabular una vez sobre [T_w,in, T_rocío + margen] y extraer [T₉, T₈] por interpolación. Elimina el perfil del lazo frío | Bajo–medio |
| **Paso adaptativo** | Paso de 2 K fijo; en líquido subenfriado (casi lineal) se puede usar 5–10 K y refinar solo entre T_burbuja y T_rocío | Medio |
| **Calcular T_burbuja y T_rocío una vez** | Delimitan exactamente los quiebres donde puede estar el mínimo de ΔT; con ellos como nodos obligatorios bastan pocos nodos intermedios | Bajo |
| **Nodos saltados** | Hoy un fallo del flash se ignora; usar la interpolación con `bubbleP/dewP` (estables en x) para rellenarlo evita huecos en la zona crítica del perfil | Medio |
| **Solo en la última pasada** | Activar el limitador únicamente cuando el lazo frío ya está cerca de converger (residuo < 10·tol) y comprobar después; las pasadas iniciales no lo necesitan con precisión | Medio (requiere validar convergencia) |
