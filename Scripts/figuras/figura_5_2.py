#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figura_4_2.py — Figura 5-2 del capítulo 5.

Desviación de las cinco curvas calculadas frente a las dos series publicadas
del trabajo de referencia, **con el interpolante analítico del apartado 4.4.1**
y no con la interpolación lineal que emplea el analizador.

El interpolante es el descrito en `interpolante_exp1.py`: la solución de
Bateman del apartado 5.1 usada como función de forma, reescalada al nivel de
los nodos de ACAB con un único factor multiplicativo k.

Uso:
    python3 figura_4_2.py CARPETA_EXP1 SERIE_COMP.csv SERIE_EXP.csv [SALIDA]

donde CARPETA_EXP1 es el directorio que contiene las cinco carpetas de caso.
"""

import glob
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LN2 = math.log(2.0)
BARN = 1.0e-24
RHO = 0.12317          # g/cm3 — masa del blanco por cm3, del fort.6

# Los cinco casos, en el orden en que se citan en el capítulo.
# La clave se busca sobre el nombre de la carpeta normalizado (sin espacios,
# puntos ni guiones y en minúsculas), de modo que funciona tanto con
# "Exp1 - v.2 - DECAY (2007)" como con "Exp1_-_v_2_-_DECAY__2007_".
CASOS = [
    ("v.1  — Haffner",              "v1",  ["v1b"], "#1b1b1b", "-"),
    ("v.1b — Haffner, flujo térm.", "v1b", [],      "#1b1b1b", "--"),
    ("v.2  — DECAY (2007)",         "v2",  [],      "#c1440e", "-"),
    ("v.3  — DECAY (2025)",         "v3",  [],      "#2e6f9e", "-."),
    ("v.4  — DECAY-XSEC (2025)",    "v4",  [],      "#4a7c3f", "-"),
]


def _normaliza(nombre):
    return "".join(c for c in nombre.lower() if c.isalnum())


def localizar(raiz, clave, excluir):
    """Carpeta de simulación cuyo nombre normalizado contiene `clave`."""
    for ruta in sorted(glob.glob(os.path.join(raiz, "*"))):
        if not os.path.isdir(ruta):
            continue
        norma = _normaliza(os.path.basename(ruta))
        if clave not in norma or any(x in norma for x in excluir):
            continue
        if os.path.exists(os.path.join(ruta, "inp.5")):
            return ruta
        hijos = [d for d in sorted(glob.glob(ruta + "/*"))
                 if os.path.exists(os.path.join(d, "inp.5"))]
        if hijos:
            return hijos[0]
    return None



# --------------------------------------------------------------------------
# Lectura de datos nucleares (los mismos ficheros que lee ACAB)
# --------------------------------------------------------------------------

def _seccion(fo, zaid, mt):
    txt = open(glob.glob(fo + "/XSECTION.dat")[0], "rb").read().decode("latin-1")
    lineas = txt.replace("\r", "").split("\n")
    for i, l in enumerate(lineas):
        if l[:8].strip() == str(zaid) and l[8:13].strip() == str(mt):
            for k in range(i + 1, i + 6):
                try:
                    return float(lineas[k].strip())
                except ValueError:
                    continue
    return None


def _decay(fo, zaid):
    txt = open(glob.glob(fo + "/DECAY.dat")[0], "rb").read().decode("latin-1")
    lineas = txt.replace("\r", "").split("\n")
    for i, l in enumerate(lineas):
        if len(l) > 12 and l[6:12].strip() == str(zaid):
            c = l[12:].split()
            return float(c[1]), float(c[6]), float(lineas[i + 1].split()[5]) / 100.0
    return None, None, None


def _nodos_acab(fo):
    """Instantes y actividad del I-131 de la tabla de enfriamiento del fort.6.

    La tabla buena es la que sigue al encabezado
    `NUCLIDE RADIOACTIVITY, DISINTEGRATIONS/SEC` **después** del bloque
    `CONCENTRATIONS AFTER IRRADIATION`: el fichero contiene varias tablas con
    columnas `INITIAL` y hay que anclarse a ese encabezado y no al primero que
    aparezca.

    Las columnas son INITIAL, SHUTDOWN y los instantes de enfriamiento; la
    primera se descarta y el apagado se toma como t = 0.
    """
    txt = open(glob.glob(fo + "/fort.6")[0], "rb").read().decode("latin-1")
    lineas = txt.replace("\r", "").split("\n")

    # La tabla buena es la de actividad (RADIOACTIVITY) del bloque posterior
    # a la irradiación (AFTER IRRADIATION). El fichero contiene seis tablas con
    # columnas INITIAL/SHUTDOWN y solo una cumple las dos condiciones.
    for i, l in enumerate(lineas):
        if not (l.strip().startswith("INITIAL") and "SHUTDOWN" in l):
            continue
        contexto = "\n".join(lineas[max(0, i - 12):i])
        if "RADIOACTIVITY" not in contexto or "AFTER IRRADIATION" not in contexto:
            continue
        ts = [0.0] + [float(x) for x in l.split()[2:]]
        for j in range(i + 1, i + 4000):
            if lineas[j][:7].strip() == "I131":
                v = [float(x) for x in lineas[j][7:].split()]
                return ts, v[1:1 + len(ts)]
    raise ValueError("no se encuentra la tabla de actividad de enfriamiento")


# --------------------------------------------------------------------------
# El interpolante analítico
# --------------------------------------------------------------------------

def interpolante(fo, tmax_h=6.0):
    """Devuelve f(t_h) -> Bq/cm3 y el factor de reescalado k.

    PASO 1  forma: la cadena de Bateman integrada con los datos de esta
            simulación, con paso de un segundo.
    PASO 2  nivel: k = media de A_ACAB / A_analitica sobre los nodos de
            enfriamiento, excluido el instante del apagado.
    PASO 3  evaluación: f(t) = k * A_analitica(t).
    """
    sg = _seccion(fo, 521300, 1020)
    sm = _seccion(fo, 521300, 1021)

    lineas = open(fo + "/inp.5", encoding="latin-1").read().replace("\r", "").split("\n")
    phi = float(lineas[[i for i, x in enumerate(lineas) if "Block #3" in x][0] + 1])
    i5 = [i for i, x in enumerate(lineas) if "Block #5" in x][0]
    n_te = float(lineas[i5 + 2].split()[0]) * 1e24

    t_g, _, _ = _decay(fo, 521310)
    t_m, fit, _ = _decay(fo, 521311)
    t_i, _, _ = _decay(fo, 531310)
    _, _, abund = _decay(fo, 521300)

    n_130 = n_te * abund
    lg, lm, li = LN2 / t_g, LN2 / t_m, LN2 / t_i
    src_g, src_m = n_130 * sg * BARN * phi, n_130 * sm * BARN * phi

    def rhs(t, y, irradiando):
        te_g, te_m, i131 = y
        a, b = (src_g, src_m) if irradiando else (0.0, 0.0)
        return [a - lg * te_g + fit * lm * te_m,
                b - lm * te_m,
                lg * te_g + (1 - fit) * lm * te_m - li * i131]

    # PASO 1
    y_eoi = solve_ivp(rhs, (0, 10.0008), [0, 0, 0], args=(True,),
                      rtol=1e-13, atol=1e-10).y[:, -1]
    malla = np.arange(0.0, tmax_h * 3600 + 1.0, 1.0)
    sol = solve_ivp(rhs, (0, malla[-1]), y_eoi, args=(False,), t_eval=malla,
                    rtol=1e-13, atol=1e-10)
    act = li * sol.y[2]
    horas = malla / 3600.0

    # PASO 2 — factor de escala, no ajuste de forma
    ts, acab = _nodos_acab(fo)
    coc = [a / np.interp(t, horas, act)
           for t, a in zip(ts, acab) if t > 0 and a > 0]
    k = float(np.mean(coc))

    # PASO 3
    return (lambda t: k * np.interp(t, horas, act)), k


# --------------------------------------------------------------------------

def leer_serie(ruta):
    """Serie digitalizada de la suite: '#' de cabecera, 't;A', coma decimal."""
    t, a = [], []
    for linea in open(ruta, encoding="utf-8-sig"):
        linea = linea.strip()
        if not linea or linea.startswith("#") or linea.lower().startswith("t;"):
            continue
        campos = linea.replace(",", ".").split(";")
        if len(campos) >= 2:
            t.append(float(campos[0]))
            a.append(float(campos[1]))
    orden = np.argsort(t)
    return np.array(t)[orden], np.array(a)[orden]


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    raiz, ruta_comp, ruta_exp = sys.argv[1], sys.argv[2], sys.argv[3]
    salida = sys.argv[4] if len(sys.argv) > 4 else "figura_4_2"

    t_comp, a_comp = leer_serie(ruta_comp)
    t_exp, a_exp = leer_serie(ruta_exp)

    fig, ejes = plt.subplots(2, 1, figsize=(7.0, 7.8), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.115, right=0.985, top=0.955, bottom=0.175)

    paneles = [
        (ejes[0], t_comp, a_comp, "frente a la serie computacional publicada"),
        (ejes[1], t_exp,  a_exp,  "frente a la serie experimental publicada"),
    ]

    resumen = {}
    etiquetas_finales = [[], []]      # (y, texto, color) por panel
    for etiqueta, clave, excluir, color, trazo in CASOS:
        fo = localizar(raiz, clave, excluir)
        if fo is None:
            print("  aviso: no se encuentra el caso %s" % etiqueta)
            continue
        f, k = interpolante(fo)
        resumen[etiqueta] = {"k": k}
        for n_panel, (eje, ts, act, _) in enumerate(paneles):
            desv = 100.0 * (f(ts) / 1e6 / RHO / act - 1.0)
            eje.plot(ts, desv, trazo, color=color, linewidth=1.5,
                     marker="o", markersize=3.2, label=etiqueta)
            resumen[etiqueta].setdefault("medias", []).append(desv.mean())
            etiquetas_finales[n_panel].append((desv[-1], desv.mean(), color))

    for n_panel, (eje, ts, _, titulo) in enumerate(paneles):
        eje.axhline(0.0, color="#999999", linewidth=0.8, zorder=0)
        eje.set_ylabel("Desviación  [%]")
        eje.grid(True, alpha=0.25, linewidth=0.6)
        eje.set_title(titulo, fontsize=10, loc="left", pad=6)
        eje.margins(x=0.02)

        # Sesgo medio de cada serie, junto a su extremo derecho. Es la columna
        # de la tabla 4-2, para que figura y tabla se lean juntas.
        x0, x1 = eje.get_xlim()
        eje.set_xlim(x0, x1 + 0.155 * (x1 - x0))
        vistos = []
        for y_fin, media, color in sorted(etiquetas_finales[n_panel]):
            while any(abs(y_fin - v) < 0.045 * (eje.get_ylim()[1] - eje.get_ylim()[0])
                      for v in vistos):
                y_fin += 0.045 * (eje.get_ylim()[1] - eje.get_ylim()[0])
            vistos.append(y_fin)
            eje.annotate("%+.2f %%" % media, xy=(x1, y_fin),
                         xytext=(x1 + 0.02 * (x1 - x0), y_fin),
                         color=color, fontsize=8.2, va="center",
                         fontweight="bold")

    ejes[1].set_xlabel("Tiempo tras el fin de la irradiación  [h]")

    # Leyenda fuera de los ejes: las dos gráficas están llenas y cualquier
    # recuadro interior tapa una serie o una anotación.
    manejadores, rotulos = ejes[0].get_legend_handles_labels()
    fig.legend(manejadores, rotulos, loc="lower center", ncol=3,
               fontsize=8.4, frameon=False, bbox_to_anchor=(0.5, 0.005),
               columnspacing=1.6, handlelength=2.4)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (salida, ext), dpi=200)
    print("Figura escrita en %s.pdf y %s.png\n" % (salida, salida))

    print("%-30s %9s %12s %12s" % ("caso", "k", "vs comput.", "vs exper."))
    for etiqueta, d in resumen.items():
        print("%-30s %9.5f %+11.3f %% %+11.3f %%"
              % (etiqueta, d["k"], d["medias"][0], d["medias"][1]))


if __name__ == "__main__":
    main()
