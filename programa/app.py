"""Panel de simulacion del Ciclo Kalina KCS-11 -- interfaz sobre kalina.py.

No modifica kalina.py/incertidumbre.py: solo declara variables (FIJO/BARRIDO),
arma `par`, corre resolver()/criterios()/clasificar()/supuestos() por cada
combinacion de la rejilla, y presenta/registra el resultado.

Ejecutar con:
    streamlit run app.py
"""
import sys
import time
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import motor_ui as mu  # noqa: E402
import sugerencias  # noqa: E402

st.set_page_config(page_title="Ciclo Kalina -- Panel de simulacion", layout="wide",
                    page_icon="🌀")

st.markdown("""
<style>
[data-testid="stMetricValue"] { font-family: "Consolas", monospace; }
.bloque-var { border-left: 3px solid #3a5a7a; padding-left: 0.8rem; margin-bottom: 0.3rem; }
code { font-size: 0.85rem; }
</style>
""", unsafe_allow_html=True)

st.title("🌀 Ciclo Kalina KCS-11 -- Panel de simulacion")
st.caption(
    "Interfaz sobre `kalina.py`/`incertidumbre.py`. El objetivo del programa es "
    "diagnosticar si un conjunto de datos es valido para un ciclo Kalina, no producir "
    "un veredicto de una sola planta -- ver `proyecto-kalina-objetivo-codigo` en memoria."
)

if "df_resultados" not in st.session_state:
    st.session_state["df_resultados"] = None
if "casos_detalle" not in st.session_state:
    st.session_state["casos_detalle"] = []
if "en_barrido_ejecutado" not in st.session_state:
    st.session_state["en_barrido_ejecutado"] = []
if "ultimo_modo" not in st.session_state:
    st.session_state["ultimo_modo"] = None
if "ultimas_filas" not in st.session_state:
    st.session_state["ultimas_filas"] = None


# ---------------------------------------------------------------------------
# Sidebar: modo de analisis y opciones de ejecucion
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("Modo de analisis")
    modo = st.radio(
        "Modo", ["EXPLORACION", "SELECCION", "FRONTERA"],
        help=(
            "EXPLORACION: todos los casos con su clasificacion e indicadores, sin "
            "condicion de aplicabilidad.\n\n"
            "SELECCION: filtra VALIDO, excluye fallas de titulo de turbina (O1), ordena "
            "por potencia neta descendente. Solo procede si BARRIDO esta unicamente en "
            "P_alta/P_baja/x_b/m_b.\n\n"
            "FRONTERA: ubica el/los cruce(s) entre 'existe regimen Kalina' y 'no existe' "
            "a lo largo de un parametro de cierre. Exige que algun cierre "
            "(HRVG/regenerador/condensador) este en BARRIDO."
        ),
    )
    st.divider()
    st.header("Ejecucion")
    timeout = st.number_input(
        "Tiempo limite por caso (s)", min_value=5, max_value=1800, value=90, step=5,
        help=(
            "Los puntos en la region de 6-8 kg/s de gas equivalente pueden tardar "
            "minutos (composicion de equilibrio nueva para el proceso -- ver "
            "decisiones-codigo-kalina.md, correccion C27). Un caso que excede este "
            "tiempo se marca TIEMPO_AGOTADO, no se descarta como invalido."
        ),
    )
    calcular_condensador = st.checkbox(
        "Calcular tambien el diagnostico del condensador (C_cold minimo)",
        value=False,
        help="Una pasada extra de evaluacion de estado() por caso -- mas lento.",
    )
    guardar_auto = st.checkbox(
        "Guardar en Excel automaticamente al terminar", value=True,
        help=f"Libro acumulativo unico: {mu.RUTA_EXCEL}",
    )
    st.divider()
    st.caption(f"Registro: `{mu.RUTA_EXCEL.relative_to(mu.RUTA_PROYECTO)}`")


# ---------------------------------------------------------------------------
# Widget generico FIJO/BARRIDO para una variable
# ---------------------------------------------------------------------------
def fila_variable(key, etiqueta, unidad, defecto, minimo, maximo):
    c1, c2, c3, c4, c5 = st.columns([2.4, 1, 1.1, 1.1, 0.7])
    c1.markdown(f"<div class='bloque-var'><b>{etiqueta}</b> <code>{unidad}</code></div>",
                unsafe_allow_html=True)
    modo_v = c2.selectbox("modo", ["FIJO", "BARRIDO"], key=f"modo_{key}",
                           label_visibility="collapsed")
    if modo_v == "FIJO":
        valor = c3.number_input("valor", value=float(defecto), key=f"valor_{key}",
                                 label_visibility="collapsed", format="%.5g")
        c4.empty()
        c5.empty()
        return {"modo": "FIJO", "valor": valor}
    vmin = c3.number_input("min", value=float(minimo), key=f"min_{key}",
                            label_visibility="collapsed", format="%.5g")
    vmax = c4.number_input("max", value=float(maximo), key=f"max_{key}",
                            label_visibility="collapsed", format="%.5g")
    n = c5.number_input("n", value=5, min_value=2, max_value=50, step=1,
                         key=f"n_{key}", label_visibility="collapsed")
    return {"modo": "BARRIDO", "min": vmin, "max": vmax, "n": n}


specs = {}

st.subheader("1. Variables de diseno")
hc1, hc2, hc3, hc4, hc5 = st.columns([2.4, 1, 1.1, 1.1, 0.7])
hc2.caption("modo"); hc3.caption("valor / min"); hc4.caption("max"); hc5.caption("n")
for key, (etiqueta, unidad, defecto, minimo, maximo) in mu.VARIABLES_BASE.items():
    specs[key] = fila_variable(key, etiqueta, unidad, defecto, minimo, maximo)

st.subheader("2. Cierres de los equipos")
tipo_cierre = {}
for eq, opciones in mu.CIERRES.items():
    st.markdown(f"**{opciones['etiqueta']}**")
    formulaciones = [k for k in opciones if k != "etiqueta"]
    cformu, _ = st.columns([2, 3])
    sel = cformu.radio(
        "Formulacion", formulaciones,
        format_func=lambda f, _o=opciones: _o[f][0],
        key=f"tipo_{eq}", horizontal=True, label_visibility="collapsed",
    )
    tipo_cierre[eq] = sel
    etiqueta_v, unidad_v, defecto_v, min_v, max_v = opciones[sel]
    specs[eq] = fila_variable(f"{eq}_{sel}", etiqueta_v, unidad_v, defecto_v, min_v, max_v)

st.subheader("3. Variables opcionales (precision adicional)")
st.caption(
    "Ninguna es obligatoria. Si no se declara, el motor asume la hipotesis de su "
    "docstring y `supuestos()` reporta la incertidumbre que eso deja -- se muestra "
    "en el detalle de cada caso."
)
opcionales_declaradas = []
for key, (etiqueta, unidad, defecto, minimo, maximo) in mu.VARIABLES_OPCIONALES.items():
    declarar = st.checkbox(f"Declarar -- {etiqueta} (`{unidad}`)",
                            key=f"declarar_{key}", value=False)
    if declarar:
        specs[key] = fila_variable(key, etiqueta, unidad, defecto, minimo, maximo)
        opcionales_declaradas.append(key)

en_barrido = mu.variables_en_barrido(specs)


def validar_modo(modo_sel, en_barrido_):
    if modo_sel == "SELECCION":
        no_permitidas = [v for v in en_barrido_ if v not in mu.VARIABLES_DISENO]
        if no_permitidas:
            return False, (
                "SELECCION exige que solo P_alta/P_baja/x_b/m_b esten en BARRIDO "
                "(especificacion_codigo_v2.md #8). Tambien estan en BARRIDO: "
                + ", ".join(no_permitidas)
            )
    if modo_sel == "FRONTERA":
        cierres_en_barrido = [eq for eq in mu.CIERRES if eq in en_barrido_]
        if not cierres_en_barrido:
            return False, (
                "FRONTERA exige que el parametro de cierre de algun equipo "
                "(HRVG/regenerador/condensador) este en BARRIDO."
            )
    return True, ""


st.divider()
n_casos = 1
for v in en_barrido:
    n_casos *= max(int(specs[v].get("n", 1)), 2)
st.caption(f"Rejilla actual: **{n_casos} caso(s)** ({len(en_barrido)} variable(s) en BARRIDO).")

ejecutar = st.button("▶ Ejecutar", type="primary")

if ejecutar:
    ok, msg = validar_modo(modo, en_barrido)
    if not ok:
        st.error(msg)
    else:
        combos = mu.construir_grid(specs)
        st.caption(
            "Corriendo toda la rejilla en un proceso compartido (conserva el cache "
            "interno de kalina.py entre casos -- solo se pierde y se reinicia si un "
            "caso concreto se cuelga mas alla del tiempo limite)."
        )
        barra = st.progress(0.0)
        estado_txt = st.empty()
        pares = [mu.combo_a_par(c, tipo_cierre, opcionales_declaradas) for c in combos]
        filas = [None] * len(combos)
        detalle = [None] * len(combos)
        t_total0 = time.time()

        def _on_caso(i, resultado, dt):
            filas[i] = mu.fila_de_caso(pares[i], combos[i]["m_b"], resultado, modo, dt)
            detalle[i] = dict(par=pares[i], m_b=combos[i]["m_b"], resultado=resultado)
            barra.progress((i + 1) / len(combos))
            estado_txt.caption(
                f"Caso {i + 1}/{len(combos)} (id {i}) -- {resultado['estado']} -- "
                f"{dt:.1f} s (transcurrido: {time.time() - t_total0:.0f} s)"
            )

        mu.resolver_grid(pares, timeout=timeout, calcular_condensador=calcular_condensador,
                         on_caso=_on_caso)
        df = pd.DataFrame(filas)
        st.session_state["df_resultados"] = df
        st.session_state["casos_detalle"] = detalle
        st.session_state["ultimo_modo"] = modo
        st.session_state["en_barrido_ejecutado"] = en_barrido
        st.session_state["ultimas_filas"] = filas
        st.success(f"{len(combos)} caso(s) resueltos en {time.time() - t_total0:.1f} s.")
        if guardar_auto:
            ruta, total = mu.guardar_en_excel(filas)
            st.info(f"Guardado en `{ruta}` -- {total} fila(s) acumuladas en total.")


# ---------------------------------------------------------------------------
# Resultados
# ---------------------------------------------------------------------------
df = st.session_state["df_resultados"]
if df is not None and not df.empty:
    modo_mostrado = st.session_state["ultimo_modo"]
    st.subheader(f"Resultados -- modo {modo_mostrado}")

    COLORES = {
        "VALIDO": "#1e5f3a", "VALIDO con advertencia": "#6b5b1e",
        "DEGENERADO": "#5a4a1e", "INVIABLE": "#6b1e1e",
        "NO_CONVERGE": "#444444", "TIEMPO_AGOTADO": "#333333", "ERROR": "#6b1e1e",
    }

    def _color_etiqueta(val):
        c = COLORES.get(val)
        return f"background-color: {c}; color: white" if c else ""

    st.dataframe(
        df.style.map(_color_etiqueta, subset=["etiqueta"]),
        use_container_width=True, height=420,
    )

    if not guardar_auto and st.session_state["ultimas_filas"]:
        if st.button("💾 Guardar esta corrida en el Excel acumulativo"):
            ruta, total = mu.guardar_en_excel(st.session_state["ultimas_filas"])
            st.info(f"Guardado en `{ruta}` -- {total} fila(s) acumuladas en total.")

    if modo_mostrado == "SELECCION":
        st.markdown("#### Ordenamiento SELECCION (VALIDO, sin fallas de O1, por W_net desc)")
        sel = mu.aplicar_seleccion(df)
        if sel.empty:
            st.warning(
                "Ningun caso de la rejilla quedo VALIDO tras el filtro de SELECCION. "
                "Esto es un resultado real (ver memoria `proyecto-kalina-objetivo-codigo`): "
                "el diagnostico es el entregable, no el veredicto."
            )
        else:
            st.dataframe(sel, use_container_width=True)
            top = sel.iloc[0]
            st.success(
                f"Optimo de la rejilla: W_net = {top['Wnet_kW']:.2f} kW, "
                f"eta = {top['eta_pct']:.2f} % -- "
                f"P_alta={top['P_alta_MPa']:.3f} MPa, P_baja={top['P_baja_MPa']:.3f} MPa, "
                f"x_b={top['x_b']:.3f}, m_b={top['m_b_kg_s']:.3f} kg/s"
            )

    if modo_mostrado == "FRONTERA":
        st.markdown("#### Frontera del regimen Kalina")
        en_barrido_ej = st.session_state["en_barrido_ejecutado"]
        cierres_en_barrido = [eq for eq in mu.CIERRES if eq in en_barrido_ej]
        col_a, col_b = st.columns([1, 2])
        eq_frontera = col_a.selectbox(
            "Parametro de cierre a analizar", cierres_en_barrido,
            format_func=lambda e: mu.CIERRES[e]["etiqueta"],
        )
        columna_cierre = mu.COLUMNA_DE[eq_frontera]
        otras_barridas = [v for v in en_barrido_ej if v != eq_frontera]
        columnas_agrupacion = [mu.COLUMNA_DE[v] for v in otras_barridas]
        if columnas_agrupacion:
            col_b.caption("Se agrupa por: " + ", ".join(columnas_agrupacion))
        else:
            col_b.caption("Solo el parametro de cierre esta en BARRIDO: un unico cruce global.")
        frontera = mu.calcular_frontera(df, columna_cierre, columnas_agrupacion)
        if frontera.empty:
            st.warning("No se pudo calcular la frontera (revisa que haya casos OK en la rejilla).")
        else:
            st.dataframe(frontera, use_container_width=True)
            st.caption(
                "Un tramo `frontera_desde/hasta` es el intervalo de la rejilla donde el "
                "regimen cambia de KALINA (VALIDO/VALIDO con advertencia) a otra cosa "
                "(DEGENERADO/INVIABLE/NO_CONVERGE); no se interpola dentro del tramo -- "
                "el proyecto usa rejilla, no optimizador, por la discontinuidad del "
                "separador (especificacion_codigo_v2.md #7). 'existe_kalina_en_todo_el_rango' "
                "en False y sin cruces significa que ningun punto de esa combinacion es Kalina."
            )

    st.divider()
    st.subheader("Detalle de un caso")
    idx = st.selectbox("Elegir caso (indice de la tabla de arriba)",
                       options=list(range(len(df))),
                       format_func=lambda i: f"{i} -- {df.iloc[i]['etiqueta']}")
    caso = st.session_state["casos_detalle"][idx]
    resultado = caso["resultado"]
    if resultado["estado"] != "OK":
        st.error(f"{resultado['estado']}: {resultado.get('mensaje', '')}")
        if resultado.get("traceback"):
            with st.expander("Traceback"):
                st.code(resultado["traceback"])
    else:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Estados del ciclo (1-10)**")
            filas_est = []
            for i in range(1, 11):
                e = resultado["est"][str(i)]
                filas_est.append({
                    "estado": i, "T (K)": round(e["T"], 2), "P (MPa)": e["P"],
                    "w": round(e["w"], 4), "h (kJ/kg)": round(e["h"], 3),
                    "s (kJ/kg-K)": round(e["s"], 4), "fase": e["fase"],
                    "q": round(e.get("q", float("nan")), 4),
                })
            st.dataframe(pd.DataFrame(filas_est), use_container_width=True, hide_index=True)
        with c2:
            st.markdown("**Criterios de admisibilidad**")
            filas_crit = [{"id": cid, "cumple": ok, "severidad": sev, "detalle": txt}
                         for cid, ok, sev, txt in resultado["crit"]]
            df_crit = pd.DataFrame(filas_crit)

            def _color_ok(val):
                return "" if val else "background-color: #6b1e1e; color: white"

            st.dataframe(df_crit.style.map(_color_ok, subset=["cumple"]),
                        use_container_width=True, hide_index=True, height=320)

        st.markdown("**Supuestos e incertidumbre declarada**")
        for s in resultado["sup"]:
            marca = "✅ dado" if s["declarada"] else "⚠️ supuesto"
            st.markdown(f"- **{s['variable']}** [{marca}] = {s['valor']}  \n  → {s['efecto']}")

        sug_lista = sugerencias.sugerencias_de(resultado["violados"], resultado["etiqueta"])
        if sug_lista:
            st.markdown("**Advertencias y sugerencias de correccion**")
            for cid, texto in sug_lista:
                st.warning(f"**{cid}** -- {texto}")
else:
    st.info("Configura las variables arriba y presiona **Ejecutar** para correr la rejilla.")
