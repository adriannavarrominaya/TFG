#!/usr/bin/env python3
"""
interpolante_exp1.py — funcion de interpolacion empleada en el apartado 4.4.

Este es el interpolante con el que se calculan los sesgos de la tabla 4-2. No
es una interpolacion de la salida de ACAB: es la solucion analitica de Bateman
del apartado 4.1 usada como FUNCION DE FORMA, reescalada al nivel de los nodos
que ACAB calcula.

Uso
---
    python3 interpolante_exp1.py CARPETA SERIE.csv [SERIE2.csv ...]

    CARPETA     directorio de simulacion (inp.5, fort.6, XSECTION.dat,
                DECAY.dat)
    SERIE.csv   serie de referencia digitalizada, formato de la suite:
                cabecera con lineas '#', fila 't;A' y pares separados por ';'
                con coma decimal

Ejemplo
-------
    python3 interpolante_exp1.py "v.2 - DECAY (2007)" \\
            exp1_computacional_normalizado.csv exp1_experimental_normalizado.csv

Requiere bateman_i131.py en el mismo directorio.
"""

import os
import sys

from bateman_i131 import cargar_caso, interpolar


# ==========================================================================
#  LA FUNCION, en tres pasos
# ==========================================================================

def construir_interpolante(carpeta, tmax_h=6.0, dt=1.0):
    """Devuelve f(t_h) -> actividad [Bq/cm3], y el factor de reescalado k.

    PASO 1 — Forma.
        Se integra la cadena de Bateman con los datos de los ficheros de esta
        simulacion, con paso de un segundo. El resultado es una curva continua
        A_ana(t) que reproduce la curvatura real de la rama de crecimiento.

    PASO 2 — Nivel.
        La curva analitica y la salida de ACAB no coinciden exactamente: entre
        ambas hay un desfase constante, identificado en 4.1 como convenio de
        conversion de atomos a actividad. Se corrige con un unico factor
        multiplicativo

            k = (1/N) * sum_j [ A_ACAB(t_j) / A_ana(t_j) ]

        promediado sobre los nodos de enfriamiento. Es un factor de escala, no
        un ajuste: no toca la forma de la curva, solo su nivel.

        El nodo t = 0 se excluye del promedio. En el apagado la actividad es
        cuatro ordenes de magnitud menor que en el maximo y su cociente esta
        dominado por el redondeo a cuatro cifras; incluirlo introduciria en k
        el ruido de un solo punto mal condicionado.

    PASO 3 — Evaluacion.
        f(t) = k * A_ana(t), interpolada linealmente sobre la malla de un
        segundo. A esa resolucion la interpolacion lineal es exacta a efectos
        practicos: el error de trocear una exponencial en tramos de un segundo
        es de orden (lambda*dt)^2/8, unos 10^-13 para esta cadena.

    Por que no interpolar directamente los nodos de ACAB. La malla de salida
    tiene paso de 0,25 h y la curva es concava en toda la rama de crecimiento,
    de modo que la recta entre dos nodos pasa siempre por debajo de la curva
    real. El sesgo resultante es sistematico y vale 0,30 puntos porcentuales
    en media, y 1,15 en el primer instante comparado.
    """
    _, cadena, acab = cargar_caso(carpeta)

    # PASO 1
    ts, act = cadena.enfriar(tmax_h, dt)

    # PASO 2
    cocientes = [a_acab / interpolar(ts, act, t_h)
                 for t_h, a_acab in zip(acab["t_h"], acab["A_i131"])
                 if t_h > 0 and a_acab > 0]
    k = sum(cocientes) / len(cocientes)

    # PASO 3
    def f(t_h):
        return k * interpolar(ts, act, t_h)

    return f, k


def sesgos(carpeta, serie, tmax_h=6.0, dt=1.0):
    """Desviaciones punto a punto de la curva calculada frente a una serie.

    Devuelve (lista de (t, desviacion en %), sesgo medio, recorrido de los
    diez ultimos, pendiente ajustada).

    La desviacion se define igual que en la suite:
        (A_calculada - A_serie) / A_serie * 100
    """
    f, _ = construir_interpolante(carpeta, tmax_h, dt)
    pts = [(t, 100.0 * (f(t) / a - 1.0)) for t, a in serie]
    d = [v for _, v in pts]
    media = sum(d) / len(d)
    recorrido = max(d[1:]) - min(d[1:])          # excluye el primer instante
    n = len(pts)
    mx = sum(t for t, _ in pts) / n
    my = sum(v for _, v in pts) / n
    pend = (sum((t - mx) * (v - my) for t, v in pts)
            / sum((t - mx) ** 2 for t, _ in pts))
    return pts, media, recorrido, pend


def leer_serie(ruta):
    """Serie de referencia digitalizada: lista de pares (t en h, A en MBq/g)."""
    out = []
    with open(ruta, encoding="utf-8-sig") as fh:
        for linea in fh:
            linea = linea.strip()
            if not linea or linea.startswith("#") or linea.lower().startswith("t;"):
                continue
            campos = linea.replace(",", ".").split(";")
            if len(campos) >= 2:
                out.append((float(campos[0]), float(campos[1])))
    return sorted(out)


# ==========================================================================

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    carpeta, series = sys.argv[1], sys.argv[2:]

    _, cadena, acab = cargar_caso(carpeta)
    f, k = construir_interpolante(carpeta)

    # Densidad del blanco, para pasar de Bq/cm3 a MBq/g igual que la suite:
    # se deduce de los propios ficheros dividiendo las dos magnitudes.
    rho = _densidad(carpeta)

    print("CASO: %s" % os.path.basename(carpeta.rstrip("/")))
    print("  factor de reescalado k = %.6f" % k)
    print("  densidad del blanco    = %.7f g/cm3" % rho)

    for ruta in series:
        serie = leer_serie(ruta)
        # La serie esta en MBq/g y el interpolante en Bq/cm3: se convierte.
        serie_bq = [(t, a * 1e6 * rho) for t, a in serie]
        pts, media, rec, pend = sesgos(carpeta, serie_bq)
        print("\n  Serie: %s  (%d puntos, %.4f a %.4f h)"
              % (os.path.basename(ruta), len(serie), serie[0][0], serie[-1][0]))
        print("    sesgo medio                  %+.3f %%" % media)
        print("    recorrido (10 ultimos)       %.3f pp" % rec)
        print("    pendiente ajustada           %+.3f pp/h" % pend)
        print("    punto a punto:")
        for t, v in pts:
            print("      t=%7.4f h   %+7.3f %%" % (t, v))


def _densidad(carpeta):
    """Masa del blanco por cm3 [g/cm3], del bloque CONCENTRATIONS(GRAM).

    Es el divisor que la suite emplea al convertir Bq/cm3 a MBq/g, y no
    coincide con el que se obtiene aplicando pesos atomicos estandar: ACAB
    compone la masa con numeros masicos ponderados por las abundancias de su
    librería

    Hay que anclar la busqueda al encabezado del bloque: el fort.6 contiene
    varias filas TOTAL antes de esa, pertenecientes a bloques distintos.
    """
    import re
    from bateman_i131 import _leer
    lineas = _leer(os.path.join(carpeta, "fort.6")).split("\n")
    for i, linea in enumerate(lineas):
        if linea.startswith(" CONCENTRATIONS(GRAM) AFTER"):
            for k in range(i, min(i + 40, len(lineas))):
                m = re.match(r"^ TOTAL\s+([\d.]+E[+-]\d+)", lineas[k])
                if m:
                    return float(m.group(1))
    raise ValueError("no se ha podido leer la masa del blanco del fort.6")


if __name__ == "__main__":
    main()
