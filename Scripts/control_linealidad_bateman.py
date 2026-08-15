#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
control_linealidad_bateman.py — verificación independiente del control Sigma R_i

QUÉ HACE
--------
Recalcula, leyendo directamente los ficheros fort.6, el control de linealidad
de las ecuaciones de Bateman del análisis de contribución por cadenas:

    R_i    = A_i(IFINAL, t*) / A_ref(IFINAL, t*)
    Sigma R_i  = 1   si se han incluido todos los isótopos iniciales

donde A_i es la actividad del isótopo objetivo en la ejecución monoisotópica
del isótopo i, y A_ref la del caso de referencia completo. La igualdad se
sigue de que las ecuaciones de Bateman son lineales en las concentraciones
iniciales y de que la matriz de transición es la misma en todas las
ejecuciones.

POR QUÉ EXISTE ESTE SCRIPT
--------------------------
La suite ya calcula Sigma R_i. Este script NO lo sustituye: lo recalcula por
un camino que no comparte código con ella. En este proyecto, todos los
errores de fondo se han detectado por discrepancia entre dos implementaciones
independientes, nunca releyendo una sola.

EL DIAGNÓSTICO QUE INCLUYE, Y ES LA PARTE MÁS ÚTIL
--------------------------------------------------
Sigma R_i = 1 es condición NECESARIA pero NO SUFICIENTE. La primera ejecución
de este control dio 0,9999 y era falsa: los inp.5 monoisotópicos llevaban
INPT = 1, con lo que ACAB expandía cada identificador a elemento natural y
los R_i resultaban ser las abundancias isotópicas. La suma valía 1 porque las
abundancias suman 1.

El delator es el cociente R_i / C_i: si sale CONSTANTE entre isótopos, se
están midiendo abundancias y no contribuciones. Físicamente es imposible que
todos los isótopos rindan lo mismo por átomo, porque no todos distan el mismo
número de capturas del objetivo. El script lo comprueba siempre.

USO
---
    python control_linealidad_bateman.py <carpeta_analisis> [--iso I131]

Espera la estructura que genera la suite:
    <carpeta>/            fort.6 de la referencia, o --ref <ruta>
    <carpeta>/iso_XXXX/   una carpeta por isótopo, cada una con su fort.6

Opciones útiles:
    --ref RUTA        fort.6 del caso de referencia, si no está en la raíz
    --instante N      índice del nodo de enfriamiento a evaluar (0 = apagado)
    --pico            evalúa en el máximo de la referencia (declara el nodo)
"""

import argparse
import glob
import os
import re
import sys

# ---------------------------------------------------------------------------
# Lectura del fort.6
# ---------------------------------------------------------------------------

FIN_DE_TABLA = ("NUMBER OF ATOMS", "PHOTON RELEASE", "CONCENTRATIONS(GRAM)")


def _leer(ruta):
    with open(ruta, encoding="latin-1") as f:
        return f.read().split("\n")


def serie_actividad(ruta_fort6, isotopo):
    """Serie de actividad del isótopo durante el enfriamiento, en Bq/cm3.

    Fusiona los conjuntos temporales (TIME SET) que ACAB escribe por separado,
    uno por tarjeta de los bloques 7 y 8, y descarta las columnas que no son
    instantes:

      - INITIAL : estado de partida, no pertenece a la serie.
      - RESTART : en los bloques de CONTINUACIÓN duplica el último punto del
                  bloque anterior y se descarta. En el PRIMER bloque del
                  enfriamiento es un punto real —el fin de la irradiación— y
                  se conserva.

    Descarta también las tablas BY ZONE, que duplican las BY INTERVAL.
    """
    sim = isotopo.upper()
    m = re.match(r"^([A-Z]+)\s*(\d+M?)$", sim)
    if not m:
        raise ValueError("isótopo mal escrito: %s (ejemplo: I131, TE131M)" % isotopo)
    elem, masa = m.group(1), m.group(2)
    patron = re.compile(r"^\s{1,3}%s\s*%s\s+(.*)$" % (elem, masa))

    lineas = _leer(ruta_fort6)
    bloques, dentro, por_zona = [], False, False

    for ln in lineas:
        if "CONCENTRATIONS" in ln and "BY ZONE" in ln:
            por_zona = True
        if "CONCENTRATIONS" in ln and "BY INTERVAL" in ln:
            por_zona = False
        if "NUCLIDE RADIOACTIVITY" in ln:
            dentro = not por_zona
            continue
        if any(t in ln for t in FIN_DE_TABLA):
            dentro = False
        if dentro:
            g = patron.match(ln)
            if g:
                try:
                    bloques.append([float(x) for x in g.group(1).split()])
                except ValueError:
                    pass

    if not bloques:
        raise RuntimeError("no se encontró %s en las tablas de actividad de %s"
                           % (isotopo, ruta_fort6))

    serie = list(bloques[0][1:])          # fuera INITIAL; RESTART aquí es el apagado
    for b in bloques[1:]:
        serie.extend(b[2:])               # fuera INITIAL y RESTART de continuación
    return serie, len(bloques)


def concentraciones_iniciales(ruta_fort6):
    """Inventario inicial del eco de la entrada: nucleido -> átomos/cm3.

    Se lee de la tabla NUMBER OF ATOMS, columna INITIAL. Es la composición YA
    EXPANDIDA a isótopos, que es lo que ACAB usó de verdad.
    """
    lineas = _leer(ruta_fort6)
    dentro, conc = False, {}
    for ln in lineas:
        if "NUMBER OF ATOMS" in ln:
            dentro = True
            continue
        if dentro and ("CONCENTRATIONS(GRAM)" in ln or "RADIOACTIVITY" in ln):
            break
        if dentro:
            g = re.match(r"^\s{1,3}([A-Z]{1,2})\s*(\d+M?)\s+(\S+)\s+\S+", ln)
            if g:
                clave = g.group(1) + g.group(2)
                try:
                    v = float(g.group(3))
                except ValueError:
                    continue
                if clave not in conc:
                    conc[clave] = v
    return conc


# ---------------------------------------------------------------------------
# El control
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("carpeta", help="carpeta raíz del análisis de cadenas")
    p.add_argument("--iso", default="I131", help="isótopo objetivo (por defecto I131)")
    p.add_argument("--ref", default=None, help="fort.6 del caso de referencia")
    p.add_argument("--instante", type=int, default=None,
                   help="índice del nodo de enfriamiento (0 = fin de irradiación)")
    p.add_argument("--pico", action="store_true",
                   help="evaluar en el máximo de la referencia")
    args = p.parse_args()

    ref = args.ref or os.path.join(args.carpeta, "fort.6")
    if not os.path.isfile(ref):
        sys.exit("no encuentro el fort.6 de referencia: %s\n"
                 "indícalo con --ref" % ref)

    serie_ref, nb = serie_actividad(ref, args.iso)
    print("REFERENCIA  %s" % ref)
    print("  conjuntos temporales fusionados: %d   instantes: %d" % (nb, len(serie_ref)))

    # --- elección del instante, declarada ---------------------------------
    if args.pico:
        mx = max(serie_ref)
        empate = [i for i, v in enumerate(serie_ref) if v == mx]
        k = empate[0]
        if len(empate) > 1:
            print("  AVISO: el máximo empata en %d instantes (índices %s) a las cuatro"
                  % (len(empate), ", ".join(map(str, empate))))
            print("         cifras que imprime el fort.6. Se evalúa en el índice %d." % k)
    elif args.instante is not None:
        k = args.instante
    else:
        k = 0
    if k >= len(serie_ref):
        sys.exit("el índice %d no existe: la serie tiene %d instantes" % (k, len(serie_ref)))

    a_ref = serie_ref[k]
    print("  instante evaluado: índice %d   A_ref(%s) = %.5e Bq/cm3\n" % (k, args.iso, a_ref))
    if a_ref == 0:
        sys.exit("la actividad de referencia es cero en ese instante")

    conc = concentraciones_iniciales(ref)

    # --- las ejecuciones monoisotópicas -----------------------------------
    carpetas = sorted(glob.glob(os.path.join(args.carpeta, "iso_*")))
    if not carpetas:
        sys.exit("no encuentro carpetas iso_* en %s" % args.carpeta)

    filas, suma = [], 0.0
    for d in carpetas:
        nombre = os.path.basename(d)[4:]
        f6 = os.path.join(d, "fort.6")
        if not os.path.isfile(f6):
            print("  %-8s sin fort.6, se omite" % nombre)
            continue
        try:
            s, _ = serie_actividad(f6, args.iso)
            a_i = s[k] if k < len(s) else float("nan")
        except RuntimeError:
            a_i = 0.0                      # el objetivo no aparece: contribución nula
        r = a_i / a_ref
        suma += r
        filas.append((nombre, conc.get(nombre, float("nan")), a_i, r))

    # --- salida ------------------------------------------------------------
    print("%-9s %14s %14s %12s %14s" % ("isótopo", "C_i [át/cm3]", "A_i [Bq/cm3]", "R_i", "R_i / C_i"))
    print("-" * 70)
    cocientes = []
    for nombre, c, a_i, r in filas:
        q = r / c if c and c == c and c != 0 else float("nan")
        if q == q and r > 0:
            cocientes.append(q)
        print("%-9s %14.4e %14.4e %12.6f %14.4e" % (nombre, c, a_i, r, q))
    print("-" * 70)
    print("%-9s %14s %14s %12.6f" % ("SUMA", "", "", suma))

    # --- veredicto ---------------------------------------------------------
    print("\nCONTROL 1 — linealidad de Bateman")
    d = abs(suma - 1.0)
    if d < 5e-4:
        print("  Sigma R_i = %.6f  ->  PASA (desviación %.1e, del orden del redondeo)" % (suma, d))
    else:
        print("  Sigma R_i = %.6f  ->  desviación %.1e" % (suma, d))
        print("  Si la selección de isótopos es PARCIAL, esto es la cobertura y no un fallo.")
        print("  Si es completa, hay que investigarlo antes de publicar nada.")

    print("\nCONTROL 2 — ¿son contribuciones o abundancias?")
    if len(cocientes) < 2:
        print("  no hay bastantes isótopos con contribución no nula para juzgar")
    else:
        disp = (max(cocientes) - min(cocientes)) / (sum(cocientes) / len(cocientes))
        print("  dispersión relativa de R_i/C_i entre isótopos: %.2e" % disp)
        if disp < 0.05:
            print("  ATENCIÓN: el cociente es CONSTANTE entre isótopos.")
            print("  Eso es físicamente imposible para la producción de un nucleido que")
            print("  dista distinto número de capturas de cada uno. El síntoma clásico es")
            print("  que los inp.5 monoisotópicos lleven INPT = 1, con lo que ACAB expande")
            print("  cada identificador a ELEMENTO NATURAL y lo que se mide son abundancias.")
            print("  Compruébalo antes de dar por bueno el Sigma R_i.")
        else:
            print("  el cociente VARÍA entre isótopos, como debe: se están midiendo")
            print("  contribuciones y no abundancias.")

    print("\nRecuerda al citar: el instante de evaluación, y que las cifras de un caso")
    print("congelado de verificación no son resultados físicos.")


if __name__ == "__main__":
    main()
