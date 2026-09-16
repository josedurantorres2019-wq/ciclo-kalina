"""Capa de orquestacion entre kalina.py y la interfaz Streamlit.

No modifica kalina.py ni incertidumbre.py: solo construye `par`, corre
resolver()/criterios()/clasificar()/supuestos() sobre una rejilla de casos
(producto cartesiano de las variables declaradas BARRIDO) y empaqueta el
resultado para la tabla en pantalla y el registro en Excel.

Cada caso corre en un PROCESO propio via subprocess (worker_caso.py o
worker_lote.py), no en un hilo ni via multiprocessing.Process: kalina.py
mantiene cache a nivel de modulo (_puros, _flash, _env, _ultimaT) sin locks,
pensada para un solo hilo de ejecucion, asi que correr casos concurrentes en
el mismo proceso arriesgaria corromperlas. subprocess (en vez de
multiprocessing con metodo 'spawn') evita ademas que Windows reimporte el
script de Streamlit como si fuera el worker -- ver docstring de worker_caso.py.
El proceso aislado tambien es lo que permite cancelar de verdad un caso que se
cuelga en la region cara conocida (6-8 kg/s de gas, ver
decisiones-codigo-kalina.md C27).
"""
import itertools
import json
import os
import queue
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

RUTA_PROYECTO = Path(__file__).resolve().parent.parent
RUTA_EXCEL = RUTA_PROYECTO / "datos" / "registro_casos_kalina.xlsx"
WORKER_SCRIPT = Path(__file__).resolve().parent / "worker_caso.py"
WORKER_LOTE_SCRIPT = Path(__file__).resolve().parent / "worker_lote.py"

_ESTE_DIR = str(Path(__file__).resolve().parent)
if _ESTE_DIR not in sys.path:
    sys.path.insert(0, _ESTE_DIR)
import sugerencias  # noqa: E402  (mismo directorio)

# ---------------------------------------------------------------------------
# Definicion de variables expuestas en la interfaz
# ---------------------------------------------------------------------------
# Variables de diseno: las UNICAS admitidas en BARRIDO bajo el modo SELECCION
# (especificacion_codigo_v2.md #8: "solo procede cuando las variables en
# BARRIDO son variables de diseno -- presiones, composicion, caudal").
VARIABLES_DISENO = ["P_alta", "P_baja", "w_b", "m_b"]

# etiqueta, unidad, valor por defecto, min sugerido, max sugerido
VARIABLES_BASE = {
    "P_alta": ("Presion alta", "MPa", 3.0, 0.5, 6.0),
    "P_baja": ("Presion baja", "MPa", 0.4, 0.1, 2.0),
    "w_b": ("Fraccion basica de amoniaco (x_b)", "-", 0.5, 0.30, 0.80),
    "m_b": ("Caudal de solucion basica", "kg/s", 1.0, 0.05, 20.0),
    "T_f": ("Temperatura de la fuente", "K", 623.15, 400.0, 800.0),
    "T_amb": ("Temperatura ambiente / sumidero", "K", 300.15, 260.0, 320.0),
    "eta_t": ("Eficiencia isentropica turbina", "-", 0.85, 0.5, 1.0),
    "eta_p": ("Eficiencia isentropica bomba", "-", 0.75, 0.5, 1.0),
}

# Cada equipo de cierre admite dos formulaciones (aproximacion o efectividad);
# el usuario elige UNA para toda la corrida.
CIERRES = {
    "cierre_HRVG": {
        "etiqueta": "Cierre del HRVG (caldera de recuperacion)",
        "dT_app": ("Approach dT_app", "K", 20.0, 0.0, 60.0),
        "eps": ("Efectividad eps_HRVG", "-", 0.85, 0.5, 0.999),
    },
    "cierre_reg": {
        "etiqueta": "Cierre del regenerador",
        "dT_pp": ("Pinch dT_pp", "K", 5.0, 0.0, 40.0),
        "eps": ("Efectividad eps_reg", "-", 0.75, 0.3, 0.999),
    },
    "cierre_cond": {
        "etiqueta": "Cierre del condensador",
        "dT_app": ("Approach dT_app", "K", 5.0, 0.0, 30.0),
        "eps": ("Efectividad eps_cond", "-", 0.80, 0.3, 0.999),
    },
}

# Variables opcionales de precision: NUNCA obligatorias (memoria
# kalina-variables-opcionales-con-incertidumbre). Si no se declaran, resolver()
# asume la hipotesis de su docstring y supuestos() reporta la incertidumbre.
VARIABLES_OPCIONALES = {
    "C_g": ("Capacidad calorifica de la fuente (C_g)",
            "kW/K por kg/s de m_b", 10.0, 0.5, 100.0),
    "dT_pp_gas": ("Margen de pinzamiento del HRVG", "K", 0.0, 0.0, 20.0),
    "T_min_gas": ("Temperatura minima del gas (rocio acido)", "K", 393.15, 350.0, 450.0),
    "C_cold": ("Capacidad calorifica del agua de enfriamiento (C_cold)",
               "kW/K por kg/s de m_b", 10.0, 0.5, 100.0),
    "T_agua_in": ("Temperatura de entrada del agua", "K", 300.15, 270.0, 320.0),
    "dT_pp_cond": ("Margen de pinzamiento del condensador", "K", 0.0, 0.0, 20.0),
    "T_agua_out_max": ("Temperatura maxima de salida del agua (vertido)", "K", 313.15, 280.0, 330.0),
    "cp_agua": ("cp del agua", "kJ/kg-K", 4.18, 3.5, 4.5),
}

ES_KALINA = {"VALIDO", "VALIDO con advertencia"}


# ---------------------------------------------------------------------------
# Rejilla: de specs (FIJO/BARRIDO por variable) a lista de combinaciones
# ---------------------------------------------------------------------------
def valores_de(spec):
    """spec: {"modo": "FIJO", "valor": x} o {"modo": "BARRIDO", "min", "max", "n"}."""
    if spec["modo"] == "FIJO":
        return [float(spec["valor"])]
    n = max(int(spec["n"]), 2)
    return [float(v) for v in np.linspace(spec["min"], spec["max"], n)]


def construir_grid(specs):
    """specs: dict variable -> spec. Devuelve la lista de combinaciones (dict
    variable -> valor), producto cartesiano de las declaradas BARRIDO."""
    claves = list(specs.keys())
    listas = [valores_de(specs[c]) for c in claves]
    return [dict(zip(claves, tupla)) for tupla in itertools.product(*listas)]


def variables_en_barrido(specs):
    return [v for v, s in specs.items() if s["modo"] == "BARRIDO"]


def combo_a_par(combo, tipo_cierre, opcionales_declaradas):
    """Arma el dict `par` que espera kalina.resolver() a partir de una
    combinacion de la rejilla. m_b no entra: resolver() trabaja normalizado a
    m_b=1, la interfaz escala las salidas despues."""
    par = {
        "P_alta": combo["P_alta"], "P_baja": combo["P_baja"], "w_b": combo["w_b"],
        "T_f": combo["T_f"], "T_amb": combo["T_amb"],
        "eta_t": combo["eta_t"], "eta_p": combo["eta_p"],
        "cierre_HRVG": {tipo_cierre["cierre_HRVG"]: combo["cierre_HRVG"]},
        "cierre_reg": {tipo_cierre["cierre_reg"]: combo["cierre_reg"]},
        "cierre_cond": {tipo_cierre["cierre_cond"]: combo["cierre_cond"]},
    }
    for opc in opcionales_declaradas:
        if opc in combo:
            par[opc] = combo[opc]
    return par


def _texto_cierre(d):
    (k, v), = d.items()
    return f"{k}={v:.4g}"


def _valor_cierre(d):
    (_, v), = d.items()
    return v


# Nombre de columna en la tabla/Excel para cada variable de la rejilla (usado
# por el modo FRONTERA para agrupar por "todo lo demas" al buscar el cruce
# sobre la variable barrida elegida como parametro de cierre).
COLUMNA_DE = {
    "P_alta": "P_alta_MPa", "P_baja": "P_baja_MPa", "w_b": "x_b", "m_b": "m_b_kg_s",
    "T_f": "T_f_K", "T_amb": "T_amb_K", "eta_t": "eta_t", "eta_p": "eta_p",
    "cierre_HRVG": "cierre_HRVG_valor", "cierre_reg": "cierre_reg_valor",
    "cierre_cond": "cierre_cond_valor",
    **{v: f"{v}_declarado" for v in VARIABLES_OPCIONALES},
}


# ---------------------------------------------------------------------------
# Ejecucion de un caso en proceso aislado (subprocess), con timeout real
# ---------------------------------------------------------------------------
def resolver_caso(par, timeout=90, calcular_condensador=False):
    """Corre un caso en `python worker_caso.py` como subprocess. Devuelve dict
    con 'estado' en {"OK", "NO_CONVERGE", "ERROR", "TIEMPO_AGOTADO"}."""
    entrada = json.dumps(dict(par=par, calcular_condensador=calcular_condensador))
    try:
        proc = subprocess.run(
            [sys.executable, str(WORKER_SCRIPT)],
            input=entrada, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return dict(estado="TIEMPO_AGOTADO", mensaje=f"excedio {timeout} s sin resolver")
    if proc.returncode != 0 or not proc.stdout.strip():
        return dict(estado="ERROR",
                    mensaje=f"el proceso worker fallo (codigo {proc.returncode})",
                    traceback=proc.stderr)
    return _empaquetar_salida(json.loads(proc.stdout))


def _empaquetar_salida(salida):
    if salida["estado"] == "ok":
        return dict(estado="OK", est=salida["est"], ind=salida["ind"], crit=salida["crit"],
                    etiqueta=salida["etiqueta"], violados=salida["violados"], sup=salida["sup"])
    return dict(estado="NO_CONVERGE" if salida["estado"] == "no_converge" else "ERROR",
                mensaje=salida.get("mensaje", ""), traceback=salida.get("traceback", ""))


def resolver_grid(pares, timeout=90, calcular_condensador=False, on_caso=None,
                  n_workers=None):
    """Resuelve una lista de `par` repartida en `n_workers` procesos
    worker_lote.py en paralelo, cada uno con su propio cache interno de
    kalina.py y su cadena de h1_semilla dentro de su bloque contiguo (doc 09,
    "Paralelismo con K workers").

    La rejilla se parte en bloques CONTIGUOS para preservar la continuidad
    numerica: la semilla h1 de un caso se alimenta con la del anterior dentro
    del mismo proceso, igual que el camino serial. Los bloques cruzan solo en
    la frontera.

    Si un caso se cuelga mas de `timeout` segundos en un worker, se mata ese
    worker (los demas siguen con cache intacto), se marca ese caso
    TIEMPO_AGOTADO, y el resto de su bloque se relanza en un proceso nuevo. Un
    worker que termina inesperadamente (crash) marca ERROR en el caso en curso
    y redistribuye el resto igual.

    `n_workers=None` usa auto-deteccion: min(K, n_casos) con K = cores - 1
    (minimo 1). `n_workers=1` reproduce el camino serial anterior. Las
    rejillas de < 8 casos corren con 1 worker (un arranque frio extra pesa mas
    que el paralelismo en rejillas chicas).

    El timeout aplica POR CASO esperado dentro de cada worker, igual que
    antes. `on_caso(i, resultado, duracion_s)` se llama apenas se conoce el
    resultado del caso `i`; en paralelo el orden de llegada no es el orden de
    la lista.
    """
    n = len(pares)
    if n == 0:
        return []

    if n_workers is None:
        cores = os.cpu_count() or 1
        n_workers = max(1, min(cores - 1, n)) if n >= 8 else 1
    else:
        n_workers = max(1, min(int(n_workers), n))

    resultados = [None] * n
    pendientes = list(range(n))
    activos = {}
    cola = queue.Queue()
    _wid_seq = [0]
    completados = [0]

    # -- helpers internos ---------------------------------------------------

    def _lanzar_worker(casos):
        wid = _wid_seq[0]
        _wid_seq[0] += 1
        proc = subprocess.Popen(
            [sys.executable, str(WORKER_LOTE_SCRIPT)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1,
        )
        entrada = json.dumps([dict(par=pares[i], calcular_condensador=calcular_condensador)
                              for i in casos])
        proc.stdin.write(entrada)
        proc.stdin.close()

        stderr_buf = []

        def _leer_stderr(stream=proc.stderr, destino=stderr_buf):
            for linea in stream:
                destino.append(linea)

        h_stderr = threading.Thread(target=_leer_stderr, daemon=True)
        h_stderr.start()

        state = dict(proc=proc, casos=list(casos), j=0, t0=time.time(),
                     stderr_buf=stderr_buf, h_stderr=h_stderr, dead=False)
        activos[wid] = state

        def _leer_stdout(stream=proc.stdout, destino=cola, _wid=wid):
            for linea in stream:
                destino.put(("line", _wid, linea))
            destino.put(("fin", _wid, None))

        threading.Thread(target=_leer_stdout, daemon=True).start()
        return wid

    def _matar_worker(wid):
        w = activos[wid]
        w["proc"].kill()
        w["proc"].wait()
        w["h_stderr"].join(timeout=2)
        w["dead"] = True

    def _marcar_tiempo_agotado(wid):
        w = activos[wid]
        if w["j"] >= len(w["casos"]):
            _matar_worker(wid)
            return
        idx = w["casos"][w["j"]]
        resultados[idx] = dict(estado="TIEMPO_AGOTADO",
                               mensaje=f"excedio {timeout} s sin resolver")
        completados[0] += 1
        if on_caso:
            on_caso(idx, resultados[idx], timeout)
        leftover = w["casos"][w["j"] + 1:]
        _matar_worker(wid)
        if leftover:
            pendientes[0:0] = leftover

    def _marcar_error(wid):
        w = activos[wid]
        if w["j"] < len(w["casos"]):
            idx = w["casos"][w["j"]]
            w["h_stderr"].join(timeout=2)
            err = "".join(w["stderr_buf"])
            resultados[idx] = dict(
                estado="ERROR",
                mensaje="el proceso worker termino sin resolver este caso",
                traceback=err,
            )
            completados[0] += 1
            if on_caso:
                on_caso(idx, resultados[idx], time.time() - w["t0"])
            leftover = w["casos"][w["j"] + 1:]
            if leftover:
                pendientes[0:0] = leftover
        _matar_worker(wid)

    def _despachar():
        while pendientes:
            alive = sum(1 for w in activos.values() if not w["dead"])
            if alive >= n_workers:
                break
            slots = n_workers - alive
            chunk_size = max(1, -(-len(pendientes) // slots))
            chunk = pendientes[:chunk_size]
            del pendientes[:chunk_size]
            _lanzar_worker(chunk)

    # -- bucle principal -----------------------------------------------------

    _despachar()

    while completados[0] < n:
        ahora = time.time()
        deadlines = [w["t0"] + timeout for w in activos.values() if not w["dead"]]
        if deadlines:
            wait = max(0.0, min(deadlines) - ahora)
        else:
            wait = 0.0

        try:
            tipo, wid, payload = cola.get(timeout=wait)
        except queue.Empty:
            ahora2 = time.time()
            for wid, w in list(activos.items()):
                if not w["dead"] and ahora2 - w["t0"] >= timeout:
                    _marcar_tiempo_agotado(wid)
            _despachar()
            continue

        if wid not in activos or activos[wid]["dead"]:
            continue

        w = activos[wid]

        if tipo == "fin":
            _marcar_error(wid)
            _despachar()
        else:
            idx = w["casos"][w["j"]]
            resultados[idx] = _empaquetar_salida(json.loads(payload))
            completados[0] += 1
            if on_caso:
                on_caso(idx, resultados[idx], time.time() - w["t0"])
            w["j"] += 1
            if w["j"] < len(w["casos"]):
                w["t0"] = time.time()
            else:
                w["proc"].wait(timeout=5)
                w["dead"] = True
                _despachar()

    return resultados


# ---------------------------------------------------------------------------
# Empaquetado de una fila de tabla/Excel
# ---------------------------------------------------------------------------
def fila_de_caso(par, m_b, resultado, modo, duracion_s):
    fila = {
        "fecha_hora": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "modo": modo,
        "P_alta_MPa": par.get("P_alta"),
        "P_baja_MPa": par.get("P_baja"),
        "x_b": par.get("w_b"),
        "m_b_kg_s": m_b,
        "T_f_K": par.get("T_f"),
        "T_amb_K": par.get("T_amb"),
        "eta_t": par.get("eta_t"),
        "eta_p": par.get("eta_p"),
        "cierre_HRVG": _texto_cierre(par["cierre_HRVG"]),
        "cierre_reg": _texto_cierre(par["cierre_reg"]),
        "cierre_cond": _texto_cierre(par["cierre_cond"]),
        "cierre_HRVG_valor": _valor_cierre(par["cierre_HRVG"]),
        "cierre_reg_valor": _valor_cierre(par["cierre_reg"]),
        "cierre_cond_valor": _valor_cierre(par["cierre_cond"]),
        "duracion_s": round(duracion_s, 2),
        "estado_ejecucion": resultado["estado"],
    }
    for opc in VARIABLES_OPCIONALES:
        fila[f"{opc}_declarado"] = par.get(opc)
    if resultado["estado"] == "OK":
        ind = resultado["ind"]
        sug = sugerencias.sugerencias_de(resultado["violados"], resultado["etiqueta"])
        fila.update({
            "etiqueta": resultado["etiqueta"],
            "criterios_violados": ",".join(resultado["violados"]),
            "degenerado": ind["degenerado"],
            "eta_pct": round(100 * ind["eta"], 4),
            "Wnet_kJ_kg": round(ind["Wnet"], 4),
            "Wnet_kW": round(ind["Wnet"] * m_b, 3),
            "Qi_kJ_kg": round(ind["Qi"], 4),
            "Qi_kW": round(ind["Qi"] * m_b, 3),
            "Qout_kW": round(ind["Qout"] * m_b, 3),
            "q2": round(ind["q2"], 4),
            "q4": round(ind["q4"], 4),
            "q9": round(ind["q9"], 4),
            "T2_K": round(ind["T2"], 2),
            "T9_K": round(ind["T9"], 2),
            "advertencias": " | ".join(t for _, t in sug),
        })
    else:
        fila.update({
            "etiqueta": resultado["estado"],
            "criterios_violados": "",
            "degenerado": None, "eta_pct": None, "Wnet_kJ_kg": None, "Wnet_kW": None,
            "Qi_kJ_kg": None, "Qi_kW": None, "Qout_kW": None,
            "q2": None, "q4": None, "q9": None, "T2_K": None, "T9_K": None,
            "advertencias": resultado.get("mensaje", "") + " | " +
                            sugerencias.ETIQUETA_SUGERENCIA.get(resultado["estado"], ""),
        })
    return fila


# ---------------------------------------------------------------------------
# Registro acumulativo en Excel, directo al vault (sin MCP de Obsidian)
# ---------------------------------------------------------------------------
def guardar_en_excel(filas):
    """Agrega filas al libro acumulativo unico resultados/registro_casos_kalina.xlsx.
    Decision del usuario (2026-09-13): un solo libro para todo el proyecto, no
    uno por sesion; no se versiona en git (ver plan-interfaz-web-kalina.md)."""
    RUTA_EXCEL.parent.mkdir(parents=True, exist_ok=True)
    nuevo = pd.DataFrame(filas)
    if RUTA_EXCEL.exists():
        existente = pd.read_excel(RUTA_EXCEL)
        combinado = pd.concat([existente, nuevo], ignore_index=True)
    else:
        combinado = nuevo
    combinado.to_excel(RUTA_EXCEL, index=False)
    return RUTA_EXCEL, len(combinado)


# ---------------------------------------------------------------------------
# Modo SELECCION: filtra VALIDO, excluye fallas de O1, ordena por Wnet desc
# ---------------------------------------------------------------------------
def aplicar_seleccion(df):
    valido = df[df["etiqueta"].isin(ES_KALINA)].copy()
    if "criterios_violados" in valido.columns:
        valido = valido[~valido["criterios_violados"].fillna("").str.contains(r"\bO1\b")]
    return valido.sort_values(["Wnet_kW", "eta_pct"], ascending=[False, False])


# ---------------------------------------------------------------------------
# Modo FRONTERA: a partir de la rejilla ya resuelta, ubica el/los cruce(s)
# entre "existe regimen Kalina" y "no existe", a lo largo del parametro de
# cierre barrido. No corre casos nuevos (decision ya tomada del proyecto:
# rejilla, no optimizador -- ver especificacion_codigo_v2.md #7).
# ---------------------------------------------------------------------------
def calcular_frontera(df, columna_cierre, columnas_agrupacion):
    """Devuelve un DataFrame con, por cada combinacion de columnas_agrupacion,
    el tramo [valor_no_kalina, valor_kalina] del parametro columna_cierre donde
    cambia si el punto es KALINA (VALIDO/VALIDO con advertencia) o no."""
    filas = []
    grupos = ([((), df)] if not columnas_agrupacion
              else df.groupby(columnas_agrupacion, dropna=False))
    for clave, grupo in grupos:
        g = grupo.sort_values(columna_cierre)
        es_k = g["etiqueta"].isin(ES_KALINA).to_numpy()
        vals = g[columna_cierre].to_numpy()
        cruces = []
        for i in range(len(vals) - 1):
            if es_k[i] != es_k[i + 1]:
                cruces.append((vals[i], vals[i + 1], es_k[i + 1]))
        base = dict(zip(columnas_agrupacion, clave if isinstance(clave, tuple) else (clave,)))
        if not cruces:
            base.update({
                "existe_kalina_en_todo_el_rango": bool(es_k.all()),
                "existe_kalina_en_algun_punto": bool(es_k.any()),
                f"{columna_cierre}_frontera_desde": None,
                f"{columna_cierre}_frontera_hasta": None,
                "lado_kalina": None,
            })
            filas.append(base)
        else:
            for lo, hi, kalina_en_hi in cruces:
                fila = dict(base)
                fila.update({
                    "existe_kalina_en_todo_el_rango": bool(es_k.all()),
                    "existe_kalina_en_algun_punto": bool(es_k.any()),
                    f"{columna_cierre}_frontera_desde": lo,
                    f"{columna_cierre}_frontera_hasta": hi,
                    "lado_kalina": "valores mayores" if kalina_en_hi else "valores menores",
                })
                filas.append(fila)
    return pd.DataFrame(filas)