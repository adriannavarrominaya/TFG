#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figura_4_3.py — Figura 4-3 del capítulo 4.

Actividad específica del yodo durante el enfriamiento, una curva por caso, con
**las dos** líneas de techo sin portador y la marca vertical en el nodo de
3,75 h.

Por qué no sale del analizador
------------------------------
La gráfica de A_esp(t) del analizador no dibuja el techo: el manual lo menciona
solo como propiedad física —«≈4,60×10⁹ MBq/g»— y no como elemento trazado. No
es una limitación, sino una decisión razonable de la herramienta, porque **el
techo no es único**: se calcula con la semivida del ¹³¹I de la biblioteca de
cada simulación, y la variante que replica los datos de la referencia sustituye
ese campo (693 400 s frente a 693 200 s).

    techo = λ·N_A / M(¹³¹I)      con M = 130,906126 u

    T½ = 693 200 s  →  4,6000×10⁹ MBq/g   (v.2, v.3, v.4)
    T½ = 693 400 s  →  4,5987×10⁹ MBq/g   (v.1, v.1b)

El dato por simulación sí existe: la exportación `I131_aesp_yodo_*.csv` lo
publica en su cabecera, una línea por caso.

Cómo se reconstruye A_esp
-------------------------
A_esp = A(¹³¹I) / masa de yodo, con la masa compuesta por **todos** los isótopos
de yodo, incluidos los estables ¹²⁷I e ¹²⁹I, que son el portador.

El ¹²⁷I es estable y su actividad es cero, de modo que su inventario **no puede
reconstruirse desde N = A/λ**: se obtiene por conservación del número másico
A = 127, restando al inventario inicial lo que aún queda como ¹²⁷Te y ¹²⁷ᵐTe.
Confundir ese punto fue el origen de una discrepancia del 1,1 % que costó tres
rondas de auditoría en el experimento 4.

Uso:
    python3 figura_4_3.py CARPETA_EXP1 [SALIDA]
"""

import glob
import math
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NA = 6.02214076e23
LN2 = math.log(2.0)

# Masas atómicas, en u. Se usa la masa atómica y no el número másico: la
# diferencia es del 0,072 % y afecta a los valores absolutos, no a las
# fracciones del techo (defecto F19 de la herramienta, corregido).
MASA = {
    "I125": 124.904630, "I126": 125.905624, "I127": 126.904473,
    "I128": 127.905809, "I129": 128.904988, "I130": 129.906674,
    "I130M": 129.906674, "I131": 130.906126, "I132": 131.907997,
    "I132M": 131.907997, "I133": 132.907797, "I133M": 132.907797,
    "I134": 133.909758, "I134M": 133.909758, "I135": 134.910059,
}

# clave normalizada, exclusiones, semivida del I-131 [s], color, trazo
CASOS = [
    ("v.1  — Haffner",              "v1",  ["v1b"], 693400.0, "#1b1b1b", "-"),
    ("v.1b — Haffner, flujo térm.", "v1b", [],      693400.0, "#1b1b1b", "--"),
    ("v.2  — DECAY (2007)",         "v2",  [],      693200.0, "#c1440e", "-"),
    ("v.3  — DECAY (2025)",         "v3",  [],      693200.0, "#2e6f9e", "-."),
    ("v.4  — DECAY-XSEC (2025)",    "v4",  [],      693200.0, "#4a7c3f", "-"),
]


def _normaliza(nombre):
    return "".join(c for c in nombre.lower() if c.isalnum())


def localizar(raiz, clave, excluir):
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
# Lectura del fort.6
# --------------------------------------------------------------------------
# Se usa `parse6.py`, el lector con el que se han verificado las cifras del
# capítulo. Conviene no reescribirlo a la ligera: el fort.6 tiene tres trampas
# encadenadas y las tres son silenciosas.
#
#   1. Cada tabla se imprime **dos veces**, «BY ZONE» y «BY INTERVAL».
#      Mezclarlas duplica instantes.
#   2. Las tablas anchas se parten en **bloques de columnas**, cada uno con su
#      propia cabecera. Quedarse con el primero pierde la mitad de los nodos.
#   3. Dentro de un bloque hay **saltos de página**: líneas que empiezan por
#      "1". Cortar ahí deja fuera los nucleidos de la segunda mitad de la tabla,
#      que es justo donde está el yodo.

import parse6


def leer_caso(carpeta):
    """Instantes de enfriamiento [h] y A_esp del yodo [MBq/g]."""
    tablas = parse6.parse(os.path.join(carpeta, "fort.6"))
    atomos, actividad = tablas["atoms_irr"], tablas["act_cool"]

    eoi = max(t for t in atomos["I131"] if isinstance(t, float))
    yodos = [n for n in actividad if n.startswith("I") and n[1:2].isdigit()]

    # λ efectiva de cada isótopo, deducida del propio fichero: λ = A/N en el
    # apagado. Evita depender de una tabla de semividas externa.
    lam = {n: actividad[n]["SHUTDOWN"] / atomos[n][eoi]
           for n in yodos
           if atomos.get(n, {}).get(eoi, 0) > 0 and actividad[n]["SHUTDOWN"] > 0}
    lam_te = {n: actividad[n]["SHUTDOWN"] / atomos[n][eoi]
              for n in ("TE127", "TE127M")
              if atomos.get(n, {}).get(eoi, 0) > 0
              and actividad.get(n, {}).get("SHUTDOWN", 0) > 0}

    n_te127_0 = sum(atomos.get(n, {}).get(eoi, 0.0) for n in ("TE127", "TE127M"))
    n_i127_0 = atomos.get("I127", {}).get(eoi, 0.0)

    def masa_yodo(clave):
        m = sum(actividad[n][clave] / lam[n] * MASA.get(n, 131.0) / NA
                for n in lam if n != "I127")
        n_te127 = sum(actividad[n][clave] / lam_te[n] for n in lam_te)
        m += (n_i127_0 + (n_te127_0 - n_te127)) * MASA["I127"] / NA
        return m

    claves = ["SHUTDOWN"] + parse6.cool_times(actividad["I131"])
    ts = np.array([0.0] + [float(c) for c in claves[1:]])
    aesp = np.array([actividad["I131"][c] / 1e6 / masa_yodo(c) for c in claves])
    return ts, aesp


def main():
    raiz = sys.argv[1] if len(sys.argv) > 1 else "."
    salida = sys.argv[2] if len(sys.argv) > 2 else "figura_4_3"

    fig, eje = plt.subplots(figsize=(7.0, 5.2))
    fig.subplots_adjust(left=0.115, right=0.885, top=0.965, bottom=0.235)

    resumen = {}
    techos = {}
    for etiqueta, clave, excluir, t12, color, trazo in CASOS:
        fo = localizar(raiz, clave, excluir)
        if fo is None:
            print("  aviso: no se encuentra %s" % etiqueta)
            continue
        ts, aesp = leer_caso(fo)
        eje.plot(ts, aesp / 1e9, trazo, color=color, linewidth=1.5,
                 marker="o", markersize=3.0, label=etiqueta)
        techo = NA / MASA["I131"] * LN2 / t12 / 1e6
        techos.setdefault(round(techo, -3), []).append(etiqueta.split()[0])
        i375 = int(np.argmin(np.abs(ts - 3.75)))
        resumen[etiqueta] = (aesp[0], aesp[i375], aesp[0] / aesp[i375],
                             100 * aesp[0] / techo, 100 * aesp[i375] / techo)

    # Las dos líneas de techo. Se separan un 0,029 %, de modo que a la escala
    # de la gráfica son indistinguibles: se traza una y la ampliación de la
    # esquina superior muestra que son dos. Enseñarlo es más honesto que
    # afirmarlo en el pie, y es el punto del apartado.
    lista = sorted(techos.items())
    for techo, casos in lista:
        eje.axhline(techo / 1e9, color="#8a3ffc", linewidth=1.1,
                    linestyle=(0, (6, 3)), alpha=0.9, zorder=1)
    eje.annotate("techo sin portador", xy=(0.02, lista[-1][0] / 1e9),
                 xycoords=("axes fraction", "data"), ha="left", va="bottom",
                 fontsize=8.6, color="#8a3ffc")

    eje.axvline(3.75, color="#999999", linewidth=0.9, linestyle=":", zorder=0)
    eje.annotate("3,75 h", xy=(3.75, 0.965), xycoords=("data", "axes fraction"),
                 fontsize=8.2, color="#666666", ha="right", va="top",
                 rotation=90)

    eje.set_xlabel("Tiempo tras el fin de la irradiación  [h]")
    eje.set_ylabel("A$_{esp}$ del yodo  [10⁹ MBq de ¹³¹I / g de yodo]")
    eje.grid(True, alpha=0.25, linewidth=0.6)
    eje.set_ylim(0, 4.95)
    eje.margins(x=0.015)

    # Eje derecho: fracción del techo, que es la magnitud que el capítulo
    # discute. Se referencia al techo de v.2/v.3/v.4; la diferencia con el otro
    # es del 0,029 % y no se aprecia en este eje.
    techo_ref = lista[-1][0]
    eje2 = eje.twinx()
    eje2.set_ylim(0, 4.95e9 / techo_ref * 100)
    eje2.set_ylabel("fracción del techo  [%]", color="#555555")
    eje2.tick_params(axis="y", colors="#555555", labelsize=8.6)

    # Ampliación: los dos techos, separados
    lupa = eje.inset_axes([0.33, 0.105, 0.42, 0.20])
    for (techo, casos), estilo in zip(lista, [(0, (5, 2)), (0, (1.6, 1.6))]):
        lupa.axhline(techo / 1e9, color="#8a3ffc", linewidth=1.3,
                     linestyle=estilo)
        lupa.annotate("%.4f×10⁹   %s" % (techo / 1e9, " y ".join(casos)),
                      xy=(0.035, techo / 1e9), xycoords=("axes fraction", "data"),
                      xytext=(0, 5), textcoords="offset points",
                      fontsize=7.2, color="#8a3ffc", va="bottom", ha="left")
    centro = sum(k for k, _ in lista) / len(lista) / 1e9
    lupa.set_ylim(centro - 0.0016, centro + 0.0022)
    lupa.set_xlim(0, 1)
    lupa.set_xticks([])
    lupa.tick_params(axis="y", labelsize=6.4)
    lupa.set_title("ampliación: son dos techos, no uno", fontsize=7.4,
                   pad=3, color="#555555")
    for lado in lupa.spines.values():
        lado.set_color("#bbbbbb")

    manejadores, rotulos = eje.get_legend_handles_labels()
    fig.legend(manejadores, rotulos, loc="lower center", ncol=3, fontsize=8.4,
               frameon=False, bbox_to_anchor=(0.5, 0.008),
               columnspacing=1.6, handlelength=2.4)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (salida, ext), dpi=200)
    print("Figura escrita en %s.pdf y %s.png\n" % (salida, salida))

    print("%-30s %12s %12s %9s %10s %10s"
          % ("caso", "A_esp EOI", "A_esp 3,75h", "cociente", "% techo", "% techo"))
    for k, v in resumen.items():
        print("%-30s %12.4E %12.4E %9.3f %9.1f%% %9.1f%%"
              % (k, v[0], v[1], v[2], v[3], v[4]))


if __name__ == "__main__":
    main()
