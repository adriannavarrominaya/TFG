#!/usr/bin/env python3
"""
bateman_i131.py — control analitico independiente de la cadena del I-131.

Resuelve las ecuaciones de Bateman de la cadena

        Te-130 --(n,gamma)--> Te-131   --beta-->  I-131
               --(n,gamma)--> Te-131m  --IT-->    Te-131
                                       --beta-->  I-131

y compara el resultado con la salida de ACAB. TODOS los datos se leen de los
ficheros de la propia simulacion: no hay ninguna constante fisica incrustada
en el codigo salvo ln(2) y el numero de Avogadro.

    sigma(Te130->Te131), sigma(Te130->Te131m)   <- XSECTION.dat
    T(Te131), T(Te131m), T(I131), FIT, FB       <- DECAY.dat
    abundancia de Te-130                        <- DECAY.dat
    densidad atomica de Te, flujo, t_irr, malla <- inp.5
    actividades e inventario de contraste       <- fort.6

Uso basico
----------
    python3 bateman_i131.py CARPETA [CARPETA ...]

Cada CARPETA es un directorio de simulacion que contiene inp.5, fort.6,
XSECTION.dat y DECAY.dat. Con varias carpetas, ademas del control individual
se imprime la comparacion entre casos.

Opciones
--------
    --tmax H        Extension del barrido fino para localizar el maximo,
                    en horas (por defecto 6).
    --dt S          Paso del barrido fino, en segundos (por defecto 1).
    --curva H1,H2   Imprime la curva interpolada en los instantes dados,
                    en horas, separados por comas. Sirve para usar la
                    solucion como interpolante continuo (apartado 4.4.1).
    --descomponer   Con exactamente dos carpetas, separa la contribucion de
                    las secciones eficaces y la de los datos de decaimiento
                    ejecutando las dos combinaciones cruzadas.
    --csv FICHERO   Vuelca la curva calculada en CSV.

Ejemplos
--------
    # control del apartado 4.1 sobre los cuatro casos
    python3 bateman_i131.py "v.1 - Haffner" "v.2 - DECAY (2007)" \\
                            "v.3 - DECAY (2025)" "v.4 - DECAY-XSECTIONS (2025)"

    # descomposicion del apartado 4.3.1
    python3 bateman_i131.py --descomponer "v.1 - Haffner" "v.2 - DECAY (2007)"

    # interpolante para los sesgos del apartado 4.4
    python3 bateman_i131.py --curva 0.2709,0.4805,0.7204 "v.2 - DECAY (2007)"
"""

import argparse
import math
import os
import re
import sys

LN2 = math.log(2.0)
N_AVOGADRO = 6.02214076e23

# Identificadores ACAB de los nucleidos de la cadena.
ZA_TE130, ZA_TE131, ZA_TE131M, ZA_I131 = "521300", "521310", "521311", "531310"
MT_NG, MT_NG_M = "1020", "1021"      # (n,gamma) a fundamental y a metaestable


# --------------------------------------------------------------------------
# Lectura de ficheros
# --------------------------------------------------------------------------

def leer_xsection(ruta):
    """Devuelve {(ZAID, MT): sigma_efectiva en barn} de XSECTION.dat.

    Cada reaccion ocupa cuatro lineas: cabecera y tres de datos, de las que
    la sigma colapsada a un grupo es el primer valor de la cuarta.
    """
    lineas = _leer(ruta).split("\n")
    sigmas = {}
    for i, linea in enumerate(lineas[:-3]):
        m = re.match(r"^ *(\d{4,7}) +(\d{2,4}) +\d+ ", linea)
        if m:
            try:
                sigmas[(m.group(1), m.group(2))] = float(lineas[i + 3].split()[0])
            except (ValueError, IndexError):
                pass
    if not sigmas:
        raise ValueError("no se ha leido ninguna seccion eficaz de %s" % ruta)
    return sigmas


def leer_decay(ruta, zaid):
    """Semivida [s], ramas y abundancia isotopica [%] de un nucleido.

    Formato de la biblioteca de decaimiento (manual ACAB v.2008): el registro
    ocupa dos lineas; THALF y FIT estan en la primera, FB y ABUN en la segunda.
    """
    lineas = _leer(ruta).split("\n")
    for i, linea in enumerate(lineas[:-1]):
        if re.match(r"^\s*\d+\s+%s\s" % zaid, linea):
            a, b = linea.split(), lineas[i + 1].split()
            return {
                "THALF": float(a[3]),   # semivida, s
                "FIT":   float(a[8]),   # rama de transicion isomerica
                "FB":    float(b[1]),   # rama beta-
                "ABUN":  float(b[5]),   # abundancia isotopica, %
            }
    raise ValueError("no se ha encontrado el nucleido %s en %s" % (zaid, ruta))


def leer_inp5(ruta):
    """Densidad atomica de Te y densidad de flujo.

    El fichero lleva comentarios `<Block #N ...` que delimitan cada bloque, y
    ambos valores se leen relativos a esa marca en lugar de por patron suelto:
    el bloque 3 tiene el flujo en su primera linea de datos, y el bloque 5 los
    identificadores de elemento en una linea y sus densidades atomicas en la
    siguiente, en el mismo orden.
    """
    lineas = _leer(ruta).split("\n")

    def datos_tras(marca, n=2):
        """Primeras n lineas de datos que siguen a un comentario de bloque."""
        for i, linea in enumerate(lineas):
            if linea.startswith("<") and marca in linea:
                out = []
                for k in range(i + 1, len(lineas)):
                    if lineas[k].startswith("<"):
                        break
                    if lineas[k].strip():
                        out.append(lineas[k].strip())
                    if len(out) == n:
                        break
                return out
        return []

    # --- flujo: bloque 3, primer valor
    flujo = None
    d = datos_tras("Block #3", 1)
    if d:
        flujo = float(d[0].split()[0])

    # --- densidad atomica de Te: bloque 5, columna del elemento 52
    dens_te = None
    d = datos_tras("Block #5", 2)
    if len(d) == 2:
        elementos = d[0].split()
        densidades = [float(x) for x in d[1].split()]
        for elem, dens in zip(elementos, densidades):
            # Identificador de elemento: Z seguido de ceros (520000 = telurio).
            if int(elem) // 10000 == 52:
                dens_te = dens
                break

    if dens_te is None or flujo is None:
        raise ValueError(
            "no se han podido leer el flujo (bloque 3) o la densidad de Te "
            "(bloque 5) de %s" % ruta)

    return {"dens_te": dens_te, "flujo": flujo}


def leer_fort6(ruta):
    """Actividades [Bq/cm3] e inventario [atomos/cm3] de la cadena.

    Devuelve la serie completa de actividad de I-131 sobre los instantes de
    enfriamiento y el inventario de los tres nucleidos en el apagado.
    """
    texto = _leer(ruta).replace("\r", "")
    lineas = texto.split("\n")

    # --- inventario en el apagado, tabla NUMBER OF ATOMS
    atomos = {}
    for i, linea in enumerate(lineas):
        if "CONCENTRATIONS DURING IRRADIATION BY INTERVAL" in linea:
            for k in range(i, min(i + 4000, len(lineas))):
                m = re.match(r"^\s{0,3}([A-Z]{1,2}\s*\d{1,3}M?)\s+(.*)", lineas[k])
                if m:
                    try:
                        v = [float(x) for x in m.group(2).split()]
                        if v[-1] > 0:
                            atomos[re.sub(r"\s+", "", m.group(1))] = v[-1]
                    except ValueError:
                        pass
                if lineas[k].startswith("1") and k > i + 5:
                    break
            break

    # --- actividades, bloques CONCENTRATIONS AFTER IRRADIATION
    bloques, i = [], 0
    while i < len(lineas):
        if lineas[i].startswith(" CONCENTRATIONS AFTER IRRADIATION BY INTERVAL"):
            j = i
            while j < len(lineas) and "INITIAL" not in lineas[j]:
                j += 1
            cabecera = lineas[j].split()
            tiempos = [float(x) for x in cabecera
                       if re.match(r"^\d\.\d+E[+-]\d+$", x)]
            datos, k = {}, j + 1
            while k < len(lineas) and not lineas[k].startswith("1"):
                m = re.match(r"^\s{0,3}([A-Z]{1,2}\s*\d{1,3}M?)\s+(.*)", lineas[k])
                if m:
                    try:
                        datos[re.sub(r"\s+", "", m.group(1))] = \
                            [float(x) for x in m.group(2).split()]
                    except ValueError:
                        pass
                k += 1
            bloques.append((tiempos, datos))
            i = k
        else:
            i += 1

    if not bloques:
        raise ValueError("no se ha encontrado ninguna tabla de actividad en %s" % ruta)

    # El primer bloque incluye la columna INITIAL; los siguientes, INITIAL y
    # RESTART. Se descartan para dejar solo apagado + enfriamiento.
    t_h = [0.0]
    serie = [bloques[0][1]["I131"][1]]
    t_h += list(bloques[0][0])
    serie += bloques[0][1]["I131"][2:]
    for tiempos, datos in bloques[1:]:
        t_h += list(tiempos)
        serie += datos["I131"][2:]

    return {"t_h": t_h, "A_i131": serie, "atomos": atomos}


def _leer(ruta):
    with open(ruta, encoding="latin-1", errors="replace") as f:
        return f.read().replace("\r", "")


# --------------------------------------------------------------------------
# Solucion de Bateman
# --------------------------------------------------------------------------

class Cadena:
    """Cadena Te-130 -> {Te-131, Te-131m} -> I-131, integrada por Runge-Kutta 4.

    Se integra en lugar de usar la forma cerrada porque la solucion analitica
    degenera cuando dos constantes de desintegracion coinciden, y porque el
    error de un RK4 con paso de un segundo esta varios ordenes por debajo del
    suelo de impresion de ACAB (cuatro cifras).
    """

    def __init__(self, sigma_g, sigma_m, t_te131, t_te131m, t_i131,
                 fit, fb, n_te130, flujo, t_irr_s):
        self.l1 = LN2 / t_te131
        self.lm = LN2 / t_te131m
        self.l3 = LN2 / t_i131
        self.fit, self.fb = fit, fb
        self.P1 = n_te130 * sigma_g * 1e-24 * flujo     # tasa a Te-131
        self.Pm = n_te130 * sigma_m * 1e-24 * flujo     # tasa a Te-131m
        self.t_irr_s = t_irr_s

    def _deriv(self, y, p1, pm):
        n1, nm, n3 = y
        return (p1 + self.fit * self.lm * nm - self.l1 * n1,
                pm - self.lm * nm,
                self.l1 * n1 + self.fb * self.lm * nm - self.l3 * n3)

    def _rk4(self, y, h, p1, pm):
        k1 = self._deriv(y, p1, pm)
        k2 = self._deriv([y[i] + h / 2 * k1[i] for i in range(3)], p1, pm)
        k3 = self._deriv([y[i] + h / 2 * k2[i] for i in range(3)], p1, pm)
        k4 = self._deriv([y[i] + h * k3[i] for i in range(3)], p1, pm)
        return [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                for i in range(3)]

    def irradiar(self, pasos=4000):
        """Estado al final de la irradiacion."""
        y, h = [0.0, 0.0, 0.0], self.t_irr_s / pasos
        for _ in range(pasos):
            y = self._rk4(y, h, self.P1, self.Pm)
        return y

    def enfriar(self, tmax_h=6.0, dt=1.0):
        """Serie de actividad de I-131 durante el enfriamiento.

        Devuelve (tiempos en horas, actividades en Bq/cm3), con el instante 0
        correspondiente al fin de la irradiacion.
        """
        y = self.irradiar()
        ts, act, t = [0.0], [self.l3 * y[2]], 0.0
        while t < tmax_h * 3600:
            y = self._rk4(y, dt, 0.0, 0.0)
            t += dt
            ts.append(t / 3600)
            act.append(self.l3 * y[2])
        return ts, act

    def inventario_eoi(self):
        """Numero de atomos de los tres nucleidos al final de la irradiacion."""
        n1, nm, n3 = self.irradiar()
        return {"TE131": n1, "TE131M": nm, "I131": n3}


# --------------------------------------------------------------------------
# Construccion del caso a partir de una carpeta
# --------------------------------------------------------------------------

def cargar_caso(carpeta, flujo=None):
    """Lee los cuatro ficheros y devuelve datos, cadena y salida de ACAB."""
    p = lambda n: os.path.join(carpeta, n)
    sig = leer_xsection(p("XSECTION.dat"))
    te130 = leer_decay(p("DECAY.dat"), ZA_TE130)
    te131 = leer_decay(p("DECAY.dat"), ZA_TE131)
    te131m = leer_decay(p("DECAY.dat"), ZA_TE131M)
    i131 = leer_decay(p("DECAY.dat"), ZA_I131)
    inp = leer_inp5(p("inp.5"))
    acab = leer_fort6(p("fort.6"))

    if flujo is None:
        flujo = inp["flujo"]

    # N(Te-130) exacto: densidad atomica del inp.5 por abundancia de la
    # biblioteca. NO se toma del eco del fort.6, que la imprime redondeada.
    n_te130 = inp["dens_te"] * 1e24 * te130["ABUN"] / 100.0

    # Duracion de la irradiacion, deducida del propio fort.6.
    t_irr_s = _t_irr(carpeta)

    datos = {
        "sigma_g": sig[(ZA_TE130, MT_NG)],
        "sigma_m": sig[(ZA_TE130, MT_NG_M)],
        "t_te131": te131["THALF"], "t_te131m": te131m["THALF"],
        "t_i131": i131["THALF"],
        "fit": te131m["FIT"], "fb": te131m["FB"],
        "abun": te130["ABUN"], "n_te130": n_te130,
        "flujo": flujo, "t_irr_s": t_irr_s,
    }
    cadena = Cadena(datos["sigma_g"], datos["sigma_m"], datos["t_te131"],
                    datos["t_te131m"], datos["t_i131"], datos["fit"],
                    datos["fb"], n_te130, flujo, t_irr_s)
    return datos, cadena, acab


def _t_irr(carpeta):
    """Duracion de la irradiacion en segundos.

    Se lee del inp.5, que es el dato que ACAB integra realmente, y no del eco
    del fort.6, que lo imprime convertido a anios con tres cifras: para un
    pulso de diez segundos ese redondeo introduce un +0,038 %, cuatro veces
    mayor que el que arrastra el propio inp.5.

    La unidad del historial se deduce comparando ambos valores, de modo que no
    hace falta interpretar IUNIT.
    """
    lineas = _leer(os.path.join(carpeta, "inp.5")).split("\n")
    t_inp = None
    for i, linea in enumerate(lineas):
        if linea.startswith("<") and re.search(r"Blocks?\s*#7", linea):
            for k in range(i + 1, len(lineas)):
                if lineas[k].startswith("<"):
                    break
                campos = lineas[k].split()
                if campos and re.match(r"^\d\.\d+E[+-]\d+$", campos[0]):
                    v = float(campos[0])
                    if v > 0:
                        t_inp = v
                        break
            break
    if t_inp is None:
        raise ValueError("no se ha podido leer el historial temporal del inp.5")

    m = re.search(r"CONCENTRATIONS AFTER ([\d.E+-]+) YEARS IRRADIATION",
                  _leer(os.path.join(carpeta, "fort.6")))
    if not m:
        raise ValueError("no se ha podido leer la duracion en el fort.6")
    t_aprox_s = float(m.group(1)) * 365.25 * 24 * 3600

    # Unidad del historial: la que hace coincidir ambos valores.
    unidades = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0,
                "a": 365.25 * 24 * 3600}
    factor = min(unidades.values(),
                 key=lambda f: abs(math.log(t_inp * f / t_aprox_s)))
    return t_inp * factor


def interpolar(ts, act, t_h):
    """Valor de la curva en un instante arbitrario, por interpolacion lineal
    sobre la malla fina de la solucion analitica (paso de un segundo)."""
    if t_h <= ts[0]:
        return act[0]
    if t_h >= ts[-1]:
        return act[-1]
    i = int(t_h * 3600 / ((ts[1] - ts[0]) * 3600))
    i = min(max(i, 0), len(ts) - 2)
    f = (t_h - ts[i]) / (ts[i + 1] - ts[i])
    return act[i] + f * (act[i + 1] - act[i])


# --------------------------------------------------------------------------
# Informes
# --------------------------------------------------------------------------

def control(carpeta, args):
    datos, cadena, acab = cargar_caso(carpeta)
    ts, act = cadena.enfriar(args.tmax, args.dt)

    print("=" * 78)
    print("CASO: %s" % os.path.basename(carpeta.rstrip("/")))
    print("=" * 78)
    print("  Datos leidos de los ficheros:")
    print("    sigma(Te130->Te131)   = %.6E b" % datos["sigma_g"])
    print("    sigma(Te130->Te131m)  = %.6E b" % datos["sigma_m"])
    print("    fraccion isomerica    = %.4f %%"
          % (100 * datos["sigma_m"] / (datos["sigma_g"] + datos["sigma_m"])))
    print("    T(Te131) / T(Te131m) / T(I131) = %.4g / %.4g / %.4g s"
          % (datos["t_te131"], datos["t_te131m"], datos["t_i131"]))
    print("    FIT / FB              = %.3f / %.3f" % (datos["fit"], datos["fb"]))
    print("    abundancia Te-130     = %.6f %%" % datos["abun"])
    print("    N(Te-130)             = %.6E at/cm3" % datos["n_te130"])
    print("    flujo                 = %.4E n/cm2 s" % datos["flujo"])
    print("    t_irr                 = %.4f s" % datos["t_irr_s"])

    # --- contraste contra el inventario
    inv = cadena.inventario_eoi()
    print("\n  Contraste contra el inventario del fort.6 (NUMBER OF ATOMS):")
    peor_inv = 0.0
    for n in ("TE131", "TE131M", "I131"):
        if n in acab["atomos"]:
            d = 100 * (inv[n] / acab["atomos"][n] - 1)
            peor_inv = max(peor_inv, abs(d))
            print("    %-7s Bateman=%.5E  ACAB=%.5E  %+.4f %%"
                  % (n, inv[n], acab["atomos"][n], d))

    # --- contraste contra la actividad
    print("\n  Contraste contra la tabla de actividad del fort.6:")
    desv = []
    for t_h, a_acab in zip(acab["t_h"], acab["A_i131"]):
        if a_acab > 0:
            desv.append(100 * (interpolar(ts, act, t_h) / a_acab - 1))
    media = sum(desv) / len(desv)
    sigma = (sum((d - media) ** 2 for d in desv) / len(desv)) ** 0.5
    print("    n = %d instantes" % len(desv))
    print("    desviacion media    = %+.4f %%" % media)
    print("    desviacion tipica   = %.4f pp" % sigma)
    print("    desviacion maxima   = %.4f %%" % max(abs(d) for d in desv))
    print("    signo               = %d de %d por debajo de ACAB"
          % (sum(1 for d in desv if d < 0), len(desv)))

    # --- instante del maximo, que ACAB no resuelve
    i = max(range(len(act)), key=lambda k: act[k])
    print("\n  Instante del maximo (la solucion analitica no esta sujeta")
    print("  al paso de malla ni al suelo de cuatro cifras de ACAB):")
    print("    t_max = %.4f h" % ts[i])
    print("    A_max = %.4E Bq/cm3" % act[i])

    if args.curva:
        print("\n  Curva interpolada en los instantes pedidos:")
        for t_h in args.curva:
            print("    t = %8.4f h   A = %.6E Bq/cm3"
                  % (t_h, interpolar(ts, act, t_h)))

    if args.csv:
        with open(args.csv, "w", encoding="utf-8") as f:
            f.write("t_h;A_Bq_cm3\n")
            for t_h, a in zip(ts, act):
                f.write("%.6f;%.6E\n" % (t_h, a))
        print("\n  Curva volcada en %s (%d puntos)" % (args.csv, len(ts)))

    return {"carpeta": carpeta, "datos": datos, "ts": ts, "act": act,
            "max": (ts[i], act[i]), "peor_inv": peor_inv,
            "media": media, "maxdesv": max(abs(d) for d in desv)}


def descomponer(carpeta_a, carpeta_b, args):
    """Separa el efecto de las secciones eficaces del de los datos de
    decaimiento, ejecutando las dos combinaciones cruzadas que ninguna
    simulacion realiza. Es el metodo del apartado 4.3.1."""
    da, _, _ = cargar_caso(carpeta_a)
    db, _, _ = cargar_caso(carpeta_b)

    def maximo(sig_g, sig_m, d):
        c = Cadena(sig_g, sig_m, d["t_te131"], d["t_te131m"], d["t_i131"],
                   d["fit"], d["fb"], d["n_te130"], d["flujo"], d["t_irr_s"])
        _, act = c.enfriar(args.tmax, args.dt)
        return max(act)

    base = maximo(db["sigma_g"], db["sigma_m"], db)
    solo_sigma = maximo(da["sigma_g"], da["sigma_m"], db)
    solo_decay = maximo(db["sigma_g"], db["sigma_m"], da)
    completo = maximo(da["sigma_g"], da["sigma_m"], da)

    print("\n" + "=" * 78)
    print("DESCOMPOSICION  %s / %s"
          % (os.path.basename(carpeta_a.rstrip("/")),
             os.path.basename(carpeta_b.rstrip("/"))))
    print("=" * 78)
    print("  solo secciones eficaces      %.5f" % (solo_sigma / base))
    print("  solo datos de decaimiento    %.5f" % (solo_decay / base))
    print("  producto de ambos            %.5f"
          % ((solo_sigma / base) * (solo_decay / base)))
    print("  sustitucion completa         %.5f" % (completo / base))
    print("  error de factorizar          %.4f pp"
          % (100 * abs((solo_sigma / base) * (solo_decay / base) - completo / base)))


def main():
    ap = argparse.ArgumentParser(
        description="Control analitico de la cadena del I-131 (apartado 4.1).")
    ap.add_argument("carpetas", nargs="+", help="carpetas de simulacion")
    ap.add_argument("--tmax", type=float, default=6.0,
                    help="extension del enfriamiento simulado, en horas")
    ap.add_argument("--dt", type=float, default=1.0,
                    help="paso de integracion en enfriamiento, en segundos")
    ap.add_argument("--curva", type=str, default=None,
                    help="instantes en horas, separados por comas")
    ap.add_argument("--descomponer", action="store_true",
                    help="con dos carpetas, separa sigma de datos de decaimiento")
    ap.add_argument("--csv", type=str, default=None,
                    help="vuelca la curva en un CSV")
    args = ap.parse_args()
    if args.curva:
        args.curva = [float(x) for x in args.curva.split(",")]

    res = [control(c, args) for c in args.carpetas]

    if len(res) > 1:
        print("\n" + "=" * 78)
        print("RESUMEN")
        print("=" * 78)
        print("  %-34s %10s %10s %10s" % ("caso", "t_max [h]", "peor inv.", "peor act."))
        for r in res:
            print("  %-34s %10.4f %9.4f%% %9.4f%%"
                  % (os.path.basename(r["carpeta"].rstrip("/")),
                     r["max"][0], r["peor_inv"], r["maxdesv"]))
        tmax = [r["max"][0] for r in res]
        print("\n  dispersion de t_max entre casos: %.4f h (%.1f min)"
              % (max(tmax) - min(tmax), 60 * (max(tmax) - min(tmax))))

    if args.descomponer:
        if len(args.carpetas) != 2:
            sys.exit("--descomponer requiere exactamente dos carpetas")
        descomponer(args.carpetas[0], args.carpetas[1], args)


if __name__ == "__main__":
    main()
