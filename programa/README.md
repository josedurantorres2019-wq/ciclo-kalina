# Ciclo Kalina KCS-11 — Simulador y validador

Simulador termodinámico del ciclo Kalina KCS-11 con mezcla amoníaco-agua
(NH₃-H₂O) sobre las propiedades de la guía **IAPWS G4-01** (Tillner-Roth y
Friend, 1998). El programa no busca emitir el veredicto de una sola planta:
**diagnostica si un conjunto de datos es válido para un ciclo Kalina** y
etiqueta cada caso (`VALIDO`, `DEGENERADO`, `INVIABLE`, `NO_CONVERGE`...) con
los criterios que incumple y la incertidumbre que dejan los supuestos.

Incluye una interfaz web (Streamlit) para correr rejillas de casos y guardar
los resultados en un Excel acumulativo.

---

## 1. Requisitos

| Componente | Versión probada | Nota |
|---|---|---|
| Python | 3.13.7 | Windows 11 |
| numpy | 2.4.6 | |
| scipy | 1.18.1 | `brentq`, `fsolve` |
| **iapws** | **1.5.5 (exacta)** | `nh3h2o.py` parchea su código fuente; con otra versión se detiene con un error |
| pandas | 3.0.5 | tablas y Excel |
| openpyxl | 3.1.5 | escritura `.xlsx` |
| streamlit | 1.63.0 | solo para la interfaz |

## 2. Instalación

Desde la carpeta `programa/`:

```bash
python -m venv .venv
```

```bash
.venv\Scripts\activate
```

```bash
pip install -r requirements.txt
```

El entorno virtual es opcional; lo único estricto es `iapws==1.5.5`
(ver [docs/01-propiedades-equilibrio-nh3h2o.md](docs/01-propiedades-equilibrio-nh3h2o.md)).

## 3. Ejecución

### 3.1 Comprobar la instalación (≈ 7 s)

Valida el motor de propiedades contra las Tablas 6, 7 y 8 de IAPWS G4-01:

```bash
python validacion_g4_01.py
```

Resultado esperado: desviación máxima ≈ `3.85e-03 %`. El único punto que
falla (burbuja x=0.6, T=500 K) está sobre el locus crítico y queda fuera del
dominio del proyecto; el propio script lo explica.

### 3.2 Interfaz web

```bash
python -m streamlit run app.py
```

Se abre en `http://localhost:8501`. También puedes hacer doble clic en
`interfaz.lnk`, que ejecuta el mismo comando (el acceso directo apunta a
`C:\Users\Usuario\AppData\Local\Programs\Python\Python313\python.exe`; si
Python está en otra ruta, edita el acceso directo o usa el comando).

Flujo en la interfaz:

1. **Modo de análisis** (barra lateral): `EXPLORACION`, `SELECCION` o `FRONTERA`.
2. **Variables de diseño**: cada una en `FIJO` (un valor) o `BARRIDO` (mín, máx, n).
3. **Cierres de los equipos**: HRVG, regenerador y condensador, por *approach/pinch* (ΔT) o por *efectividad* (ε).
4. **Variables opcionales**: ninguna es obligatoria; si no se declaran, el motor asume una hipótesis y la reporta.
5. **Ejecutar**: resuelve el producto cartesiano de los barridos.

Los resultados se agregan a `../datos/registro_casos_kalina.xlsx`
(un único libro acumulativo para todo el proyecto).

### 3.3 Pruebas de validación del ciclo

```bash
python prueba4.py
```

Reproduce un punto publicado de Elsayed et al. (2013) (η = 11.38 %) por el
método A, con sensibilidades a las hipótesis no publicadas. Tarda ≈ 8 min
(medido: 459 s). Resultado esperado: η calculada 11.18 %, `REPRODUCE`;
η_Carnot 24.13 % y η_II 46.35 % (publicadas: 24 % y 47 %).

```bash
python prueba5.py
```

Consistencia cruzada método A ↔ método B sobre el mismo punto, más un control
que **debe fallar** (demuestra que la prueba discrimina).

```bash
python incertidumbre.py
```

Estudio 1 (capacidad calorífica finita de la fuente) y Estudio 2 (banda ±0.01
de composición de G4-01). Es el más lento: varios minutos.

### 3.4 Benchmark de rendimiento

```bash
python benchmark/bench.py
```

Mide el arranque (importaciones, parches de `iapws`, primera llamada, costo de
relanzar un worker) y, por escenario, `resolver()` en frío, en caliente y en un
punto vecino, además del reparto de tiempo por algoritmo con cProfile. Cada
medición corre en un proceso limpio. Deja `reporte.md`, los JSON y los `.prof` en
`benchmark/resultados/<fecha>/`. La corrida completa tarda del orden de una
hora. La línea base medida y la lista priorizada de optimizaciones están en
[benchmark/ANALISIS.md](benchmark/ANALISIS.md). Para acotar la corrida:

```bash
python benchmark/bench.py --escenarios A_elsayed B_kalina --sin-perfil
```

### 3.5 Uso como biblioteca

```python
import kalina as k

par = dict(P_alta=3.0, P_baja=0.4, w_b=0.50,      # MPa, MPa, fraccion MASICA
           T_f=623.15, T_amb=300.15,              # K
           eta_t=0.85, eta_p=0.75,
           cierre_HRVG=dict(eps=0.85),            # o dict(dT_app=20.0)
           cierre_reg=dict(eps=0.75),             # o dict(dT_pp=5.0)
           cierre_cond=dict(eps=0.80))            # o dict(dT_app=5.0)

est, ind = k.resolver(par)                 # 10 estados + indicadores
crit = k.criterios(est, ind, par)          # lista (id, cumple, severidad, texto)
etiqueta, violados = k.clasificar(crit)
k.imprimir_supuestos(k.supuestos(est, ind, par))
print(etiqueta, violados, ind['eta'], ind['Wnet'])
```

Convenciones: T en K, P en MPa, h en kJ/kg, s en kJ/kg·K; caudales
normalizados a `m_b = 1 kg/s` (las potencias en kJ/kg equivalen a kW por kg/s
de solución básica). En la interfaz y en `par` las composiciones son
**másicas**; el motor de propiedades trabaja en **molares** y la conversión es
explícita.

## 4. Estructura

```
programa/
├── nh3h2o.py            Propiedades IAPWS G4-01 + equilibrio líquido-vapor (capa 0)
├── kalina.py            Estados, resolvedor del ciclo, pinzamiento, criterios (capa 1)
├── incertidumbre.py     Estudios de incertidumbre sobre kalina.resolver()
├── motor_ui.py          Rejilla de casos, workers, Excel, modos SELECCION/FRONTERA
├── worker_lote.py       Proceso persistente que resuelve una lista de casos (usado)
├── worker_caso.py       Proceso que resuelve un caso aislado (alternativa)
├── sugerencias.py       Criterio violado -> sugerencia de corrección
├── app.py               Interfaz Streamlit
├── validacion_g4_01.py  Validación del motor de propiedades
├── prueba4.py           Reproducción de Elsayed et al. (2013), método A
├── prueba5.py           Consistencia cruzada A <-> B
├── interfaz.lnk         Acceso directo a la interfaz
├── requirements.txt
├── benchmark/           bench.py + resultados/<fecha>/reporte.md
└── docs/                Documentación de algoritmos (ver índice)
```

Dependencias entre módulos:

```
app.py ─► motor_ui.py ─► worker_lote.py ─► kalina.py ─► nh3h2o.py ─► iapws
                    └──► sugerencias.py
incertidumbre.py, prueba4.py, prueba5.py ─► kalina.py
```

## 5. Índice de la documentación (`docs/`)

Cada documento explica **qué hace** un algoritmo, **la física o matemática que
aplica** (con sus fórmulas), **cuánto cuesta** y **dónde se puede optimizar**.
Están ordenados de la capa más baja (propiedades) a la más alta (interfaz):
leerlos en orden sigue el camino de una llamada.

| # | Documento | Algoritmo | Física / matemática |
|---|---|---|---|
| 01 | [Propiedades y equilibrio NH₃-H₂O](docs/01-propiedades-equilibrio-nh3h2o.md) | Raíz de densidad, burbuja/rocío por sustitución sucesiva + pulido de Newton, parches a `iapws` | Energía de Helmholtz, igualdad de fugacidades, valores K |
| 02 | [Flash a (T, P) y estado de la mezcla](docs/02-flash-tp-y-estado.md) | `flash_TP` con arranque tibio, `estado` | Equilibrio de fases, regla de la palanca |
| 03 | [Inversión h→T y s→T](docs/03-inversion-estado-de.md) | `estado_de`: Brent con intervalo memorizado | Monotonía de h(T), s(T); convergencia de Brent |
| 04 | [Resolución del ciclo KCS-11](docs/04-resolucion-ciclo-kcs11.md) | Lazo exterior (Brent sobre h₁) y lazo frío (punto fijo) | Balances de masa y energía, turbina/bomba isentrópicas, separador |
| 05 | [Pinzamiento en el HRVG](docs/05-pinzamiento-hrvg.md) | `pinzamiento`: recorte por mínimo acumulado | Contracorriente con capacidad finita, 2.ª ley |
| 06 | [Condensador con agua finita](docs/06-condensador-capacidad-finita.md) | `perfil_condensador`, `pinzamiento_condensador` | Espejo del HRVG del lado frío |
| 07 | [Criterios y clasificación](docs/07-criterios-clasificacion.md) | `criterios`, `clasificar` | Carnot, 2.ª ley por equipo, balance global |
| 08 | [Dimensionamiento a la fuente](docs/08-dimensionamiento-fuente.md) | `dimensionar`: raíces en m_b | Escalado C_g/m_b, rocío ácido |
| 09 | [Rejilla y ejecución por lotes](docs/09-rejilla-y-ejecucion-lotes.md) | Producto cartesiano, worker persistente, timeout, frontera | Continuación numérica (arranque tibio) |

**Resumen de optimizaciones con mayor retorno** (detalle en cada documento):

1. Jacobiano analítico y Newton en lugar de `fsolve` con diferencias finitas en el flash (doc 01–02): es el núcleo de todas las evaluaciones.
2. Newton con c_p = ∂h/∂T como paso en `estado_de` (doc 03): reduce las evaluaciones de `estado()` por inversión.
3. Aceleración de Aitken/Wegstein en el lazo frío (doc 04).
4. Memoizar `evalua(m_b)` en `dimensionar` (doc 08): hoy resuelve el ciclo dos veces por los mismos puntos.
5. Varios workers persistentes en paralelo y orden serpenteante de la rejilla (doc 09).
