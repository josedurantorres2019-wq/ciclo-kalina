"""Mapeo criterio de admisibilidad -> sugerencia de correccion accionable.

Alcance inicial (decidido con el usuario, 2026-09-13): solo los criterios que
ya aparecieron en hallazgos reales del proyecto -- O1, O2, S9, C1/C1b/C3 --
mas O3, O5, O6, S10 y N1 por ser directamente accionables con las mismas
perillas de diseno. El resto de los 20 criterios de
`contenido/especificacion_codigo_v2.md` #10 no tiene entrada todavia: se
amplia con evidencia de uso real, no de una vez -- ver
`decisiones/plan-interfaz-web-kalina.md`.
"""

SUGERENCIAS = {
    "O1": "Titulo de vapor bajo a la salida de la turbina (riesgo de erosion de alabes). "
          "Sube T_f, baja P_alta, o relaja el cierre del HRVG (mas dT_app o mas eps_HRVG) "
          "para sobrecalentar mas el estado 2.",
    "O2": "El condensador no logra condensar del todo el estado 9 (riesgo de cavitacion de la "
          "bomba aguas abajo). Sube P_baja, baja T_amb/T_agua_in del sumidero frio, aumenta su "
          "caudal (declara o sube C_cold), o relaja el cierre del condensador (mas dT_app o mas "
          "eps_cond).",
    "O3": "Potencia neta no positiva: el ciclo consume mas de lo que genera. Revisa eta_t/eta_p, "
          "o si P_alta/P_baja dejan muy poco salto entalpico disponible en la turbina.",
    "O5": "Estado 2 fuera de la banda bifasica declarada (titulo q2 muy cerca de 0 o 1): el "
          "separador no separa de forma significativa -- ciclo DEGENERADO, no un KCS-11 genuino. "
          "Ajusta P_alta o x_b, o el cierre del HRVG, para llevar el estado 2 al centro de la "
          "campana de equilibrio.",
    "O6": "La temperatura de salida del gas cae por debajo del limite de rocio acido declarado "
          "(T_min_gas). Baja el caudal de solucion basica (m_b), sube T_min_gas si el "
          "combustible lo permite, o reduce la extraccion de calor del HRVG.",
    "S9": "El perfil de temperaturas del HRVG se cruza con el del gas en un punto interior (no "
          "solo en los extremos). Baja el caudal de gas exigido (sube C_g), reduce m_b, o relaja "
          "el cierre declarado (dT_app/eps_HRVG) hasta que el margen de pinzamiento se cumpla.",
    "S10": "El perfil del condensador se cruza con el del agua de enfriamiento en un punto "
           "interior. Sube el caudal de agua (declara o sube C_cold), baja T_agua_in, o relaja "
           "el cierre del condensador.",
    "C1": "La composicion de la corriente rica o pobre no queda ordenada alrededor de x_b como "
          "exige la separacion (w_p < w_b < w_r falla). Revisa P_alta/x_b: el punto puede estar "
          "fuera del rango donde el separador opera fisicamente.",
    "C1b": "La diferencia entre la composicion rica y la pobre es menor que la banda de "
           "incertidumbre de las propiedades (+-0.01, G4-01): la separacion es marginal y el "
           "resultado no se distingue del ruido del modelo. Alejate del punto de separacion nula "
           "(sube P_alta o ajusta x_b).",
    "C3": "Las fracciones de caudal m_r/m_p salen fuera de [0,1]: revisa que P_alta/P_baja/x_b "
          "sean fisicamente consistentes entre si.",
    "N1": "El resolvedor no encontro cambio de signo en el intervalo de busqueda de h1: revisa "
          "que P_alta > P_baja y T_f > T_amb, o amplia/desplaza el rango de barrido.",
}

ETIQUETA_SUGERENCIA = {
    "NO_CONVERGE": "El resolvedor no convergio numericamente. Revisa la consistencia fisica del "
                   "punto (P_alta > P_baja, T_f > T_amb) antes de desconfiar del codigo.",
    "TIEMPO_AGOTADO": "El caso excedio el tiempo limite configurado. Es un costo conocido en la "
                       "region de composiciones de equilibrio nuevas para el proceso (ver "
                       "decisiones-codigo-kalina.md, correccion C27) -- no implica que el punto "
                       "sea invalido, solo que tarda minutos en resolver.",
    "ERROR": "El motor lanzo una excepcion no clasificada. Revisa el mensaje y el traceback antes "
             "de descartar el punto.",
}


def sugerencias_de(violados, etiqueta):
    """Lista de (criterio_o_etiqueta, texto) para los criterios violados con sugerencia
    conocida, mas la sugerencia general de la etiqueta cuando aplica (NO_CONVERGE /
    TIEMPO_AGOTADO / ERROR)."""
    out = []
    if etiqueta in ETIQUETA_SUGERENCIA:
        out.append((etiqueta, ETIQUETA_SUGERENCIA[etiqueta]))
    for cid in violados:
        if cid in SUGERENCIAS:
            out.append((cid, SUGERENCIAS[cid]))
    return out
