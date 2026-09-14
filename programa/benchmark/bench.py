"""Benchmark del motor Kalina: tiempos de arranque y consumo por algoritmo.

Mide tres cosas, cada una en un PROCESO LIMPIO (los caches de kalina.py son
globales de modulo: medir en un proceso ya usado daria tiempos "tibios"):

  1. ARRANQUE   importaciones (numpy, scipy, iapws, parches de nh3h2o, kalina,
                pandas, streamlit), primera llamada a estado() y costo de
                relanzar un worker (lo que paga motor_ui al matar un proceso).
  2. PARED      por escenario, sin perfilador: resolver() en frio, repetido
                con caches llenos (caliente), un punto vecino con h1_semilla
                (lo que ocurre en un barrido), diagnostico, pico de memoria
                y tamano de los caches.
  3. PERFIL     por escenario, resolver() en frio bajo cProfile: llamadas,
                tiempo propio y acumulado de cada algoritmo de docs/.
                cProfile encarece las llamadas pequenas (prop, _phir): usar
                los porcentajes para ORDENAR, y los tiempos de PARED como
                valor absoluto.

Uso (desde programa/):
    python benchmark/bench.py                          # todo
    python benchmark/bench.py --escenarios A_elsayed B_kalina
    python benchmark/bench.py --sin-perfil             # solo arranque + pared
    python benchmark/bench.py --solo-reporte benchmark/resultados/<carpeta>

Salida: benchmark/resultados/<fecha>/ con reporte.md, un JSON por medicion y
un .prof por escenario (abrir con `python -m pstats` o snakeviz).
"""
import argparse
import cProfile
import json
import os
import pstats
import statistics
import subprocess
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

AQUI = Path(__file__).resolve().parent
PROGRAMA = AQUI.parent
if str(PROGRAMA) not in sys.path:
    sys.path.insert(0, str(PROGRAMA))

# ---------------------------------------------------------------------------
# ESCENARIOS
# ---------------------------------------------------------------------------
# Caso base del proyecto (incertidumbre.BASE): metodo B, rama degenerada.
BASE_B = dict(P_alta=3.0, P_baja=0.4, w_b=0.50, T_f=623.15, T_amb=300.033,
              eta_t=0.85, eta_p=0.75,
              cierre_HRVG=dict(eps=0.85), cierre_reg=dict(eps=0.75),
              cierre_cond=dict(eps=0.80))
# Punto de prueba4 (Elsayed 2013), metodo A. P_baja fijada al valor que
# prueba4 obtiene con ng.bubbleP(287 K, x(w=0.55)), para que el tiempo medido
# sea solo el de resolver().
PAR_A = dict(P_alta=1.5, P_baja=0.2735498695811723, w_b=0.55, T_f=373.0,
             T_amb=283.0, eta_t=0.80, eta_p=0.80,
             cierre_HRVG=dict(dT_app=4.0), cierre_reg=dict(dT_pp=4.0),
             cierre_cond=dict(dT_app=4.0))


def _copia(par, **cambios):
    p = {k: (dict(v) if isinstance(v, dict) else v) for k, v in par.items()}
    p.update(cambios)
    return p


ESCENARIOS = {
    "A_elsayed": dict(
        desc="Metodo A, regimen Kalina (prueba4). G(h1) constante: sin Brent exterior.",
        par=PAR_A),
    "B_kalina": dict(
        desc="Metodo B sobre el mismo punto (prueba5): lazo exterior + lazo frio en regimen Kalina.",
        par=_copia(PAR_A, cierre_HRVG=dict(eps=0.8836434587827767),
                   cierre_reg=dict(eps=0.9530663541136849),
                   cierre_cond=dict(eps=0.9545779632623427))),
    "B_base": dict(
        desc="Caso base del proyecto (3.0/0.4 MPa, w_b=0.50, 350 C), metodo B, rama degenerada.",
        par=BASE_B),
    "B_base_Cg": dict(
        desc="B_base + fuente finita C_g=11 kW/K (10 kg/s de gas), dT_pp_gas=10 K: limitador del HRVG activo.",
        par=_copia(BASE_B, C_g=11.0, dT_pp_gas=10.0)),
    "B_base_Ccold": dict(
        desc="B_base + agua finita C_cold=10 kW/K: perfil del condensador en cada pasada del lazo frio.",
        par=_copia(BASE_B, C_cold=10.0)),
}

# (etiqueta, doc, sufijo de archivo, funcion). Los metodos parcheados de iapws
# viven en archivos virtuales "<parche ...>" (ver nh3h2o._aplicar_correcciones).
ALGORITMOS = [
    ("resolver (total)", "04", "kalina.py", "resolver"),
    ("tramo(h1) lazo exterior", "04", "kalina.py", "tramo"),
    ("lazo_frio", "04", "kalina.py", "lazo_frio"),
    ("cierre_hrvg", "04/05", "kalina.py", "cierre_hrvg"),
    ("cierre_reg", "04", "kalina.py", "cierre_reg"),
    ("cierre_cond", "04/06", "kalina.py", "cierre_cond"),
    ("perfil_HRVG", "05", "kalina.py", "perfil_HRVG"),
    ("pinzamiento HRVG", "05", "kalina.py", "pinzamiento"),
    ("perfil_condensador", "06", "kalina.py", "perfil_condensador"),
    ("pinzamiento_condensador", "06", "kalina.py", "pinzamiento_condensador"),
    ("estado_de (inversion h/s->T)", "03", "kalina.py", "estado_de"),
    ("estado", "02", "kalina.py", "estado"),
    ("flash_TP", "02", "kalina.py", "flash_TP"),
    ("flash_TP.F (fugacidades)", "02", "kalina.py", "F"),
    ("_mono (monofasico)", "02", "kalina.py", "_mono"),
    ("T_sat (bubbleT/dewT)", "07", "kalina.py", "T_sat"),
    ("bubbleP", "01", "nh3h2o.py", "bubbleP"),
    ("dewP", "01", "nh3h2o.py", "dewP"),
    ("_flashP (sust. sucesiva)", "01", "nh3h2o.py", "_flashP"),
    ("_polish (Newton)", "01", "nh3h2o.py", "_polish"),
    ("_satT", "01", "nh3h2o.py", "_satT"),
    ("Tsat_pure", "01", "nh3h2o.py", "Tsat_pure"),
    ("Psat_pure", "01", "nh3h2o.py", "Psat_pure"),
    ("rho_TPx (raiz de densidad)", "01", "nh3h2o.py", "rho_TPx"),
    ("P_of (P(rho))", "01", "nh3h2o.py", "P_of"),
    ("prop (propiedades)", "01", "nh3h2o.py", "prop"),
    ("iapws _prop", "01", "ammonia.py", "_prop"),
    ("iapws _phir mezcla (parcheado)", "01", "<parche _phir>", "_phir"),
    # H2ONH3._phir construye IAPWS95() y NH3() en cada llamada y pide a cada
    # puro su Helmholtz residual con TODAS las derivadas (MEoS._phir).
    ("iapws MEoS._phir (puros H2O+NH3)", "01", "iapws95.py:1924", "_phir"),
    ("iapws MEoS.__init__ (construir IAPWS95/NH3)", "01", "iapws95.py:496", "__init__"),
    ("iapws _Dphir (parcheado)", "01", "<parche _Dphir>", "_Dphir"),
    ("scipy brentq", "-", "_zeros_py.py", "brentq"),
    ("scipy fsolve", "-", "_minpack_py.py", "fsolve"),
    ("criterios", "07", "kalina.py", "criterios"),
    ("supuestos", "07", "kalina.py", "supuestos"),
]


# ---------------------------------------------------------------------------
# UTILIDADES
# ---------------------------------------------------------------------------
def pico_memoria_mb():
    """Pico de memoria residente del proceso [MB] (Windows: PeakWorkingSetSize)."""
    try:
        import ctypes
        from ctypes import wintypes as wt

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wt.DWORD), ("PageFaultCount", wt.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        k32 = ctypes.WinDLL("kernel32")
        psapi = ctypes.WinDLL("psapi")
        k32.GetCurrentProcess.restype = wt.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(PMC), wt.DWORD]
        c = PMC()
        c.cb = ctypes.sizeof(c)
        psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb)
        return c.PeakWorkingSetSize / 2**20
    except Exception:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        except Exception:
            return None


def tamano_caches():
    import kalina as k
    import nh3h2o as ng
    return dict(
        flash_presiones=len(k._flash),
        flash_soluciones=sum(len(v) for v in k._flash.values()),
        puros=len(k._puros), env=len(k._env), ultimaT=len(k._ultimaT),
        psat=len(ng._psat_cache))


def _t():
    return time.perf_counter()


# ---------------------------------------------------------------------------
# MEDICION 1: ARRANQUE (se lanza en procesos nuevos)
# ---------------------------------------------------------------------------
_SNIP_MOTOR = r"""
import json, sys, time
sys.path.insert(0, r"{prog}")
t = [time.perf_counter()]
import numpy;                          t.append(time.perf_counter())
import scipy.optimize;                 t.append(time.perf_counter())
import iapws.ammonia, iapws.iapws95;   t.append(time.perf_counter())
import nh3h2o;                         t.append(time.perf_counter())
import kalina;                         t.append(time.perf_counter())
e = kalina.estado(400.0, 3.0, 0.5);    t.append(time.perf_counter())
e = kalina.estado(400.0, 3.0, 0.5);    t.append(time.perf_counter())
e = kalina.estado(401.0, 3.0, 0.5);    t.append(time.perf_counter())
kalina.estado_de(3.0, 0.5, h=e['h'] + 20.0); t.append(time.perf_counter())
kalina.estado_de(3.0, 0.5, h=e['h'] + 21.0); t.append(time.perf_counter())
n = ["import numpy", "import scipy.optimize", "import iapws", "import nh3h2o (parches)",
     "import kalina", "estado() 1a llamada (frio)", "estado() mismo punto",
     "estado() T vecina (+1 K)", "estado_de() 1a inversion", "estado_de() inversion vecina"]
print(json.dumps({{a: t[i + 1] - t[i] for i, a in enumerate(n)}}))
"""

_SNIP_UI = r"""
import json, sys, time
sys.path.insert(0, r"{prog}")
t = [time.perf_counter()]
import pandas;     t.append(time.perf_counter())
import streamlit;  t.append(time.perf_counter())
import motor_ui;   t.append(time.perf_counter())
n = ["import pandas", "import streamlit", "import motor_ui"]
print(json.dumps({{a: t[i + 1] - t[i] for i, a in enumerate(n)}}))
"""


def medir_arranque(reps=3):
    py = sys.executable
    res = dict(reps=reps, proceso={}, motor={}, ui={})
    procesos = {
        "python vacio": ["-c", "pass"],
        "python + import kalina (relanzar worker)": ["-c", f"import sys; sys.path.insert(0, r'{PROGRAMA}'); import kalina"],
        "python + import motor_ui/streamlit": ["-c", f"import sys; sys.path.insert(0, r'{PROGRAMA}'); import motor_ui, streamlit"],
    }
    for nombre, args in procesos.items():
        ts = []
        for _ in range(reps):
            t0 = _t()
            subprocess.run([py, *args], check=True, capture_output=True)
            ts.append(_t() - t0)
        res["proceso"][nombre] = dict(mediana=statistics.median(ts), min=min(ts), max=max(ts))
        print(f"  {nombre:<45} {statistics.median(ts):7.3f} s", flush=True)
    for clave, snip in (("motor", _SNIP_MOTOR), ("ui", _SNIP_UI)):
        muestras = []
        for _ in range(reps):
            out = subprocess.run([py, "-c", snip.format(prog=PROGRAMA)], check=True,
                                 capture_output=True, text=True).stdout
            muestras.append(json.loads(out.strip().splitlines()[-1]))
        for paso in muestras[0]:
            ts = [m[paso] for m in muestras]
            res[clave][paso] = dict(mediana=statistics.median(ts), min=min(ts), max=max(ts))
            print(f"  {paso:<45} {statistics.median(ts):7.3f} s", flush=True)
    return res


# ---------------------------------------------------------------------------
# MEDICIONES 2 y 3: WORKER (un escenario en un proceso limpio)
# ---------------------------------------------------------------------------
def worker_pared(nombre):
    t0 = _t()
    import kalina as k
    t_imp = _t() - t0
    par = ESCENARIOS[nombre]["par"]
    r = dict(escenario=nombre, import_s=t_imp, pico_mem_import_mb=pico_memoria_mb())
    try:
        t0 = _t()
        est, ind = k.resolver(_copia(par))
        r["frio_s"] = _t() - t0
        r["resultado"] = dict(eta=float(ind["eta"]), degenerado=bool(ind["degenerado"]),
                              iter_frio=int(ind["iter_frio"]), h1=float(ind["h1"]))
        t0 = _t()
        crit = k.criterios(est, ind, par)
        r["etiqueta"] = k.clasificar(crit)[0]
        k.supuestos(est, ind, par, calcular_condensador=False)
        r["diagnostico_s"] = _t() - t0
        t0 = _t()
        k.diagnostico_condensador(est, ind, par)
        r["diagnostico_condensador_s"] = _t() - t0
        r["caches_tras_frio"] = tamano_caches()
        r["pico_mem_frio_mb"] = pico_memoria_mb()
        t0 = _t()
        k.resolver(_copia(par))
        r["caliente_s"] = _t() - t0
        vec = _copia(par, T_f=par["T_f"] + 2.0, h1_semilla=float(ind["h1"]))
        t0 = _t()
        try:
            k.resolver(vec)
            r["vecino_s"] = _t() - t0
        except Exception as ex:
            r["vecino_s"] = _t() - t0
            r["vecino_error"] = f"{type(ex).__name__}: {ex}"
        r["caches_final"] = tamano_caches()
    except Exception as ex:
        r["error"] = f"{type(ex).__name__}: {ex}"
        r["traceback"] = traceback.format_exc()
    r["pico_mem_final_mb"] = pico_memoria_mb()
    return r


def _es_archivo(ruta, linea, sufijo):
    """sufijo: nombre de archivo, '<parche ...>' o 'archivo.py:linea' cuando el
    mismo nombre de funcion aparece dos veces en el archivo."""
    if sufijo.startswith("<"):
        return ruta == sufijo
    arch, _, lin = sufijo.partition(":")
    return os.path.basename(ruta) == arch and (not lin or int(lin) == linea)


def worker_perfil(nombre, ruta_prof):
    import kalina as k
    par = ESCENARIOS[nombre]["par"]
    prof = cProfile.Profile()
    r = dict(escenario=nombre)
    t0 = _t()
    prof.enable()
    try:
        k.resolver(_copia(par))
    except Exception as ex:
        r["error"] = f"{type(ex).__name__}: {ex}"
        r["traceback"] = traceback.format_exc()
    finally:
        prof.disable()
    r["pared_con_perfil_s"] = _t() - t0
    prof.dump_stats(str(ruta_prof))
    st = pstats.Stats(prof)
    total = st.total_tt
    r["total_perfil_s"] = total
    filas = []
    for (f, linea, fn), (cc, nc, tt, ct, _) in st.stats.items():
        filas.append(dict(archivo=f if f.startswith("<") else os.path.basename(f),
                          ruta=f, linea=linea, funcion=fn, llamadas=nc,
                          propio_s=tt, acumulado_s=ct))
    algo = []
    for etiqueta, doc, arch, fn in ALGORITMOS:
        m = [x for x in filas if x["funcion"] == fn and _es_archivo(x["ruta"], x["linea"], arch)]
        if not m:
            continue
        ll = sum(x["llamadas"] for x in m)
        ac = sum(x["acumulado_s"] for x in m)
        pr = sum(x["propio_s"] for x in m)
        algo.append(dict(algoritmo=etiqueta, doc=doc, llamadas=ll, acumulado_s=ac,
                         propio_s=pr, pct_acumulado=100 * ac / total if total else 0,
                         pct_propio=100 * pr / total if total else 0,
                         ms_por_llamada=1000 * ac / ll if ll else 0))
    r["algoritmos"] = algo
    for x in filas:
        x.pop("ruta")
    r["top_propio"] = sorted(filas, key=lambda x: -x["propio_s"])[:30]
    return r


# ---------------------------------------------------------------------------
# ORQUESTADOR
# ---------------------------------------------------------------------------
def lanzar(modo, nombre, carpeta, timeout):
    salida = carpeta / f"{nombre}_{modo}.json"
    cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", modo,
           "--escenario", nombre, "--out", str(salida)]
    t0 = _t()
    with open(carpeta / f"{nombre}_{modo}.log", "w", encoding="utf-8") as log:
        try:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
                           cwd=str(PROGRAMA))
        except subprocess.TimeoutExpired:
            json.dump(dict(escenario=nombre, error=f"TIEMPO_AGOTADO ({timeout} s)"),
                      open(salida, "w", encoding="utf-8"))
    if not salida.exists():
        json.dump(dict(escenario=nombre, error="el worker termino sin escribir resultado (ver .log)"),
                  open(salida, "w", encoding="utf-8"))
    print(f"  [{modo:6}] {nombre:<14} {_t() - t0:8.1f} s", flush=True)


def _leer(carpeta, nombre):
    p = carpeta / nombre
    return json.load(open(p, encoding="utf-8")) if p.exists() else None


def _f(v, fmt="{:.2f}"):
    return "-" if v is None else fmt.format(v)


def escribir_reporte(carpeta):
    L = []
    meta = _leer(carpeta, "meta.json") or {}
    L.append("# Benchmark del motor Kalina\n")
    L.append(f"Fecha: {meta.get('fecha', '-')} · Python {meta.get('python', '-')} · "
             f"{meta.get('plataforma', '-')}\n")
    L.append("Generado por `benchmark/bench.py`. Cada medicion corre en un proceso "
             "limpio. Los tiempos de **pared** son los absolutos; los del **perfil** "
             "(cProfile) estan inflados por el instrumentado y sirven para ordenar.\n")

    arr = _leer(carpeta, "arranque.json")
    if arr:
        L.append("## 1. Arranque\n")
        L.append(f"Mediana de {arr['reps']} repeticiones, cada una en un proceso nuevo.\n")
        L.append("| Paso | Mediana [s] | Min [s] | Max [s] |\n|---|---:|---:|---:|")
        for grupo in ("proceso", "motor", "ui"):
            for paso, v in arr[grupo].items():
                L.append(f"| {paso} | {v['mediana']:.3f} | {v['min']:.3f} | {v['max']:.3f} |")
        L.append("")

    nombres = meta.get("escenarios", list(ESCENARIOS))
    L.append("## 2. Tiempos de pared por escenario\n")
    L.append("- **frio**: primer `resolver()` del proceso (caches vacios).\n"
             "- **caliente**: el mismo caso otra vez (caches llenos).\n"
             "- **vecino**: T_f + 2 K con `h1_semilla` (situacion de un barrido).\n")
    L.append("| Escenario | Etiqueta | Frio [s] | Caliente [s] | Vecino [s] | Frio/caliente | "
             "Diagnostico [s] | Diag. condensador [s] | Pico memoria [MB] | Soluciones flash en cache |")
    L.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for n in nombres:
        p = _leer(carpeta, f"{n}_pared.json")
        if not p:
            continue
        if "error" in p and "frio_s" not in p:
            L.append(f"| {n} | **{p['error']}** | | | | | | | | |")
            continue
        razon = (p["frio_s"] / p["caliente_s"]) if p.get("caliente_s") else None
        L.append(f"| {n} | {p.get('etiqueta', '-')} | {_f(p.get('frio_s'), '{:.1f}')} | "
                 f"{_f(p.get('caliente_s'), '{:.1f}')} | {_f(p.get('vecino_s'), '{:.1f}')} | "
                 f"{_f(razon, '{:.1f}x')} | {_f(p.get('diagnostico_s'), '{:.2f}')} | "
                 f"{_f(p.get('diagnostico_condensador_s'), '{:.1f}')} | "
                 f"{_f(p.get('pico_mem_final_mb'), '{:.0f}')} | "
                 f"{p.get('caches_final', {}).get('flash_soluciones', '-')} |")
    L.append("")
    for n in nombres:
        p = _leer(carpeta, f"{n}_pared.json")
        if p and (p.get("vecino_error") or p.get("error")):
            L.append(f"- {n}: {p.get('error') or 'vecino -> ' + p['vecino_error']}")
    L.append("")

    perfiles = {n: _leer(carpeta, f"{n}_perfil.json") for n in nombres}
    perfiles = {n: v for n, v in perfiles.items() if v and v.get("algoritmos")}
    if perfiles:
        L.append("## 3. Reparto por algoritmo (cProfile, resolver() en frio)\n")
        L.append("% del tiempo total perfilado. **Acumulado** incluye lo que llama la "
                 "funcion; **propio** es solo su cuerpo (lo que se ahorra optimizando "
                 "esa funcion sin tocar las que llama).\n")
        cab = "| Algoritmo | Doc | " + " | ".join(f"{n} acum. %" for n in perfiles) + \
              " | " + " | ".join(f"{n} propio %" for n in perfiles) + " |"
        L.append(cab)
        L.append("|---|---|" + "---:|" * (2 * len(perfiles)))
        for etiqueta, doc, _, _ in ALGORITMOS:
            acum, prop_ = [], []
            for n, v in perfiles.items():
                a = next((x for x in v["algoritmos"] if x["algoritmo"] == etiqueta), None)
                acum.append(_f(a and a["pct_acumulado"], "{:.1f}"))
                prop_.append(_f(a and a["pct_propio"], "{:.1f}"))
            if all(x == "-" for x in acum):
                continue
            L.append(f"| {etiqueta} | {doc} | " + " | ".join(acum) + " | " + " | ".join(prop_) + " |")
        L.append("")

        for n, v in perfiles.items():
            L.append(f"### {n}\n")
            L.append(f"{ESCENARIOS[n]['desc']}  \nPared con perfil: {v['pared_con_perfil_s']:.1f} s · "
                     f"total perfilado: {v['total_perfil_s']:.1f} s"
                     + (f" · **{v['error']}**" if v.get("error") else "") + "\n")
            L.append("| Algoritmo | Llamadas | Acumulado [s] | % | Propio [s] | % | ms/llamada |")
            L.append("|---|---:|---:|---:|---:|---:|---:|")
            for a in sorted(v["algoritmos"], key=lambda x: -x["acumulado_s"]):
                L.append(f"| {a['algoritmo']} | {a['llamadas']:,} | {a['acumulado_s']:.2f} | "
                         f"{a['pct_acumulado']:.1f} | {a['propio_s']:.2f} | {a['pct_propio']:.1f} | "
                         f"{a['ms_por_llamada']:.3f} |")
            L.append("\nTop 15 por tiempo propio (todas las funciones, incluidas numpy/scipy/iapws):\n")
            L.append("| Funcion | Archivo:linea | Llamadas | Propio [s] | % |")
            L.append("|---|---|---:|---:|---:|")
            for x in v["top_propio"][:15]:
                L.append(f"| `{x['funcion']}` | {x['archivo']}:{x['linea']} | {x['llamadas']:,} | "
                         f"{x['propio_s']:.2f} | {100 * x['propio_s'] / v['total_perfil_s']:.1f} |")
            L.append("")
        for n, v in {n: _leer(carpeta, f"{n}_perfil.json") for n in nombres}.items():
            if v and v.get("error") and not v.get("algoritmos"):
                L.append(f"- {n} (perfil): {v['error']}")
    (carpeta / "reporte.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return carpeta / "reporte.md"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenarios", nargs="+", choices=list(ESCENARIOS), default=list(ESCENARIOS))
    ap.add_argument("--sin-perfil", action="store_true")
    ap.add_argument("--sin-arranque", action="store_true")
    ap.add_argument("--timeout", type=float, default=3600, help="limite por proceso [s]")
    ap.add_argument("--solo-reporte", type=Path)
    ap.add_argument("--carpeta", type=Path,
                    help="agregar escenarios a una carpeta de resultados existente")
    ap.add_argument("--worker", choices=["pared", "perfil"], help=argparse.SUPPRESS)
    ap.add_argument("--escenario", help=argparse.SUPPRESS)
    ap.add_argument("--out", type=Path, help=argparse.SUPPRESS)
    a = ap.parse_args()

    if a.worker:
        os.chdir(PROGRAMA)
        if a.worker == "pared":
            r = worker_pared(a.escenario)
        else:
            r = worker_perfil(a.escenario, a.out.with_suffix(".prof"))
        a.out.write_text(json.dumps(r, indent=1), encoding="utf-8")
        return

    if a.solo_reporte:
        print(escribir_reporte(a.solo_reporte))
        return

    import platform
    if a.carpeta:
        carpeta = a.carpeta
        meta = _leer(carpeta, "meta.json") or {}
        previos = meta.get("escenarios", [])
        meta["escenarios"] = previos + [n for n in a.escenarios if n not in previos]
    else:
        carpeta = AQUI / "resultados" / datetime.now().strftime("%Y%m%d_%H%M%S")
        carpeta.mkdir(parents=True)
        meta = dict(fecha=datetime.now().isoformat(timespec="seconds"),
                    python=platform.python_version(), plataforma=platform.platform(),
                    escenarios=a.escenarios)
    json.dump(meta, open(carpeta / "meta.json", "w"))
    print(f"Resultados en {carpeta}", flush=True)
    if not a.sin_arranque:
        print("ARRANQUE", flush=True)
        json.dump(medir_arranque(), open(carpeta / "arranque.json", "w"), indent=1)
        escribir_reporte(carpeta)
    # Secuencial a proposito: correr procesos en paralelo contamina los tiempos.
    for n in a.escenarios:
        print(f"ESCENARIO {n}", flush=True)
        lanzar("pared", n, carpeta, a.timeout)
        escribir_reporte(carpeta)
        if not a.sin_perfil:
            lanzar("perfil", n, carpeta, a.timeout)
            escribir_reporte(carpeta)
    print(f"Reporte: {escribir_reporte(carpeta)}", flush=True)


if __name__ == "__main__":
    main()
