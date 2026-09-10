#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROPAGACIÓN ANALÍTICA DEL ENFRIAMIENTO — bloque 5.3 del Conjunto 2
==================================================================

Qué hace y por qué existe
-------------------------
El cálculo paramétrico de historial temporal NO se ejecuta. Los 79 nodos de irradiación de la
simulación de referencia SON el cálculo paramétrico: cada uno es el inventario que quedaría al
apagar tras irradiar ese tiempo, y desde cualquiera de ellos el enfriamiento se
propaga sin volver a ejecutar ACAB.

La razón es que sin flujo la ecuación de inventario es homogénea:

    dN/dt = A·N ,   A[i,j] = f_(j→i)·λ_j  ,  A[j,j] = −λ_j

con coeficientes constantes. La solución es exacta, N(t) = exp(A·t)·N₀, y no hay
integración numérica que controlar. Setenta y nueve historiales de irradiación
distintos, con su enfriamiento completo cada uno, salen de una sola ejecución.

Alcance del sistema
-------------------
46 especies: Te, I y Xe de A = 127 a 136 con sus isómeros, las presentes en el
DECAY.dat. Cinco ramas salen del conjunto (Xe→Cs) y NINGUNA entra: se ha comprobado
que el inventario del apagado no contiene Sb ni Sn de A ≥ 127 en cantidad que opere.
Esa es la condición que hace exacta la propagación — no que no haya salidas.

Entradas
--------
Un directorio con la simulación de referencia, que debe contener:
    inp.5        malla temporal (bloques #7/#8) y composición
    fort.6       tablas de inventario (irradiación) y de actividad (enfriamiento)
    DECAY.dat    semividas y fracciones de rama

Salidas
-------
1. VERIFICACIÓN contra los 120 nodos de enfriamiento que ACAB sí imprime.
2. CONTROL DE RAMA: qué pasa si se omite la β⁻ del ¹³¹ᵐTe (firma de rama ausente).
3. TABLA DEL CÁLCULO PARAMÉTRICO: t_máx, t_cruce, separación, precio de la espera y purezas
   para cada uno de los 79 nodos.
4. CONTROL ASINTÓTICO (E3-bis): A_esp debe decaer como el ¹³¹I cuando la masa de
   yodo se congela. Refuta el propagador, no la física.

Convenios declarados
--------------------
* MASA ATÓMICA y no número másico, en el techo y en el denominador de A_esp
  (defecto F19 del Conjunto 1). Con números másicos el techo sale 4,596704×10⁹ en
  lugar de 4,600000×10⁹ y toda la columna de A_esp queda baja un 0,073 %. Las
  purezas y fracciones no dependen del convenio.
* Las columnas del fort.6 se emparejan con su tiempo POR ORDEN y no por etiqueta:
  el fichero rotula a 3 cifras y produce colisiones.
* A partir del 2.º conjunto temporal, ACAB repite el último nodo del anterior en una
  columna RESTART. No es un nodo nuevo y se descuenta.
* El fort.6 imprime 4 cifras significativas. El suelo relativo va del 0,005 % (mantisa
  9,999) al 0,050 % (mantisa 1,000), y una comparación entre dos impresos hereda la
  suma de los dos.

Uso
---
    python3 propagacion_enfriamiento.py <directorio_de_referencia> [--csv salida.csv]

Requiere numpy y scipy.
"""

import argparse
import math
import os
import re
import sys

import numpy as np
from scipy.linalg import expm

# ---------------------------------------------------------------------------
# Constantes y datos fijos
# ---------------------------------------------------------------------------

N_AVOGADRO = 6.02214076e23

# Unidades de THALF del DECAY.dat (campo IU), en segundos. IU = 6 es estable.
UNIDADES_THALF = {1: 1.0, 2: 60.0, 3: 3600.0, 4: 86400.0, 5: 3.1557600e7,
                  7: 3.1557600e10, 8: 3.1557600e13, 9: 3.1557600e16}

SIMBOLO = {51: 'SB', 52: 'TE', 53: 'I', 54: 'XE', 55: 'CS'}

# Masas atómicas (AME2020) de los isótopos de yodo, en u. Se usan en el
# denominador de A_esp y en el techo sin portador. Ver "Convenios declarados".
MASA_ATOMICA_YODO = {
    'I127': 126.904472, 'I128': 127.905809, 'I129': 128.904984,
    'I130': 129.906671, 'I131': 130.906126, 'I132': 131.907994,
    'I133': 132.907797, 'I134': 133.909759, 'I135': 134.910048,
    'I136': 135.914604,
}

UMBRAL_PUREZA = 0.999          # umbral de calidad farmacéutica para el t_cruce


# ---------------------------------------------------------------------------
# 1. Lectura del DECAY.dat
# ---------------------------------------------------------------------------

def _campo(cadena):
    cadena = cadena.strip()
    return float(cadena) if cadena else 0.0


def leer_decay(ruta):
    """Lee el DECAY.dat con los anchos de campo del manual.

    Línea #1: I4, I8, I3, 4X, 6E10.3 -> NLB NUCL IU THALF FBX FPEC FPECX FA FIT
    Línea #2: I4, 5X, 4E10.3, ...    -> NLB FB FSF FN ...

    NUCL = Z*10^4 + A*10 + M
    """
    libreria = {}
    lineas = open(ruta, encoding='latin-1').read().replace('\r\n', '\n').split('\n')
    i = 1                                        # se salta la tarjeta de título
    while i + 1 < len(lineas):
        l1, l2 = lineas[i], lineas[i + 1]
        if not l1.strip() or l1.strip() == '-1':
            break
        nucl = int(l1[4:12])
        iu = int(l1[12:15])
        c = 19                                   # I4 + I8 + I3 + 4X
        thalf, fbx, fpec, fpecx, fa, fit = [
            _campo(l1[c + 10 * k: c + 10 * (k + 1)]) for k in range(6)]
        c2 = 9                                   # I4 + 5X
        fb, fsf, fn = [_campo(l2[c2 + 10 * k: c2 + 10 * (k + 1)]) for k in range(3)]

        z, resto = divmod(nucl, 10000)
        a, m = divmod(resto, 10)
        estable = (iu == 6)
        t_s = math.inf if estable else thalf * UNIDADES_THALF.get(iu, 1.0)
        lam = 0.0 if (estable or t_s in (0.0, math.inf)) else math.log(2) / t_s

        libreria[nucl] = dict(Z=z, A=a, M=m, iu=iu, thalf_s=t_s, lam=lam,
                                FB=fb, FBX=fbx, FPEC=fpec, FPECX=fpecx,
                                FA=fa, FIT=fit, FSF=fsf, FN=fn)
        i += 2
    return libreria


def nombre(nucl):
    z, resto = divmod(nucl, 10000)
    a, m = divmod(resto, 10)
    return SIMBOLO.get(z, 'Z%d' % z) + str(a) + ('M' * m if m else '')


def hijos(nucl, r):
    """Ramas de decaimiento -> [(nucleido_hijo, fracción)].

    Convenio del manual: FB, FPEC, FA, FIT, FSF y FN son fracciones del total y
    suman 1. FBX y FPECX NO son ramas adicionales: son la parte de FB y de FPEC
    que deja al hijo en su primer estado isomérico.

    Omitir ese reparto es un error real y documentado: deja el ¹³¹ᵐXe un 65 % bajo.
    """
    z, a = r['Z'], r['A']
    salida = []
    if r['FB'] > 0:
        fx = r['FBX']
        salida.append(((z + 1) * 10000 + a * 10, r['FB'] * (1 - fx)))
        if fx > 0:
            salida.append(((z + 1) * 10000 + a * 10 + 1, r['FB'] * fx))
    if r['FPEC'] > 0:
        fx = r['FPECX']
        salida.append(((z - 1) * 10000 + a * 10, r['FPEC'] * (1 - fx)))
        if fx > 0:
            salida.append(((z - 1) * 10000 + a * 10 + 1, r['FPEC'] * fx))
    if r['FA'] > 0:
        salida.append(((z - 2) * 10000 + (a - 4) * 10, r['FA']))
    if r['FIT'] > 0:
        salida.append((z * 10000 + a * 10, r['FIT']))          # isómero -> fundamental
    if r['FSF'] > 0:
        salida.append((z * 10000 + a * 10 + 1, r['FSF']))
    if r['FN'] > 0:
        salida.append(((z + 1) * 10000 + (a - 1) * 10, r['FN']))
    return [(h, f) for h, f in salida if f > 0]


# ---------------------------------------------------------------------------
# 2. Lectura del inp.5 y del fort.6
# ---------------------------------------------------------------------------

def leer_malla(ruta_inp5):
    """Devuelve (tiempos_irradiación, tiempos_enfriamiento) en horas, en orden."""
    lineas = open(ruta_inp5, encoding='latin-1').read().replace('\r\n', '\n').split('\n')
    i = 0
    while 'Blocks #7' not in lineas[i]:
        i += 1
    i += 1
    irradiacion, enfriamiento = [], []
    while i < len(lineas) and not lineas[i].startswith('<Block #9'):
        linea = lineas[i]
        if linea.startswith('<'):
            i += 1
            continue
        cabecera = linea.split()
        n_irr, n_total = int(cabecera[0]), int(cabecera[1])
        valores = []
        while len(valores) < n_total:
            i += 1
            valores += [float(x) for x in lineas[i].split()]
        (irradiacion if n_irr > 0 else enfriamiento).extend(valores)
        i += 1
    return irradiacion, enfriamiento


def _bloques(ruta_fort6):
    """Itera (titulo, ambito, nombres_de_columna, {nucleido: [valores]})."""
    ambito = None
    with open(ruta_fort6, encoding='latin-1', errors='replace') as fh:
        for linea in fh:
            l = linea.rstrip('\r\n')
            if 'BY ZONE' in l or 'BY INTERVAL' in l:
                ambito = 'ZONE' if 'BY ZONE' in l else 'INTERVAL'
                continue
            if 'NUMBER OF ATOMS' in l or 'NUCLIDE RADIOACTIVITY' in l:
                titulo = 'ATOMS' if 'NUMBER OF ATOMS' in l else 'ACT'
                columnas, datos = None, {}
                for l2 in fh:
                    l2 = l2.rstrip('\r\n')
                    if columnas is None:
                        if l2.strip().startswith('INITIAL'):
                            columnas = l2.split()
                        continue
                    if not l2.strip() or l2[:1] in '01':
                        break
                    clave = l2[:7].strip()
                    try:
                        datos[clave] = [float(x) for x in l2[7:].split()]
                    except ValueError:
                        break
                if columnas:
                    yield titulo, ambito, columnas, datos


def _tabla(ruta_fort6, tiempos, titulo, especies):
    """Ensambla una tabla del fort.6 emparejando columnas con tiempos por orden."""
    salida_t, series = [0.0], {}
    k = 0
    for tit, ambito, columnas, datos in _bloques(ruta_fort6):
        if tit != titulo or ambito != 'ZONE':
            continue
        hay_restart = 'RESTART' in columnas
        n = len(columnas) - (2 if hay_restart else 1)
        desplazamiento = 2 if hay_restart else 1
        bloque_t = tiempos[k:k + n]
        k += n
        for clave, valores in datos.items():
            if especies and clave not in especies:
                continue
            if clave not in series:
                # instante inicial: INITIAL en irradiación, RESTART en enfriamiento
                series[clave] = [valores[1] if hay_restart and titulo == 'ACT'
                                 else valores[0]]
            series[clave].extend(valores[desplazamiento:desplazamiento + n])
        salida_t.extend(bloque_t)
    return salida_t, series


def tabla_atomos(ruta_fort6, ruta_inp5, especies=None):
    """Inventario de átomos en cada nodo de irradiación (el 0 va delante)."""
    irr, _ = leer_malla(ruta_inp5)
    return _tabla(ruta_fort6, irr, 'ATOMS', especies)


def tabla_actividad(ruta_fort6, ruta_inp5, especies=None):
    """Actividad (des/s) en cada nodo de enfriamiento; el instante 0 es el apagado."""
    _, enf = leer_malla(ruta_inp5)
    return _tabla(ruta_fort6, enf, 'ACT', especies)


# ---------------------------------------------------------------------------
# 3. El sistema de decaimiento
# ---------------------------------------------------------------------------

def especies_TeIXe(libreria, amin=127, amax=136):
    """Te, I y Xe en el rango de masas, con isómeros, presentes en la libreria."""
    seleccion = [n for n, r in libreria.items()
                 if r['Z'] in (52, 53, 54) and amin <= r['A'] <= amax]
    return sorted(seleccion, key=lambda n: (libreria[n]['A'], libreria[n]['Z'],
                                            -libreria[n]['M']))


def matriz_decaimiento(libreria, especies):
    """A[i,j] = producción de i por decaimiento de j; A[j,j] = −λ_j.

    Las ramas que salen del conjunto se pierden (sumidero), que es correcto: lo que
    haría inexacta la propagación serían entradas, no salidas.
    """
    indice = {n: i for i, n in enumerate(especies)}
    A = np.zeros((len(especies), len(especies)))
    fugas = {}
    for n in especies:
        r = libreria[n]
        j = indice[n]
        A[j, j] = -r['lam']
        perdida = 0.0
        for h, f in hijos(n, r):
            if h in indice:
                A[indice[h], j] += f * r['lam']
            else:
                perdida += f
        if perdida > 1e-9:
            fugas[nombre(n)] = (perdida, r['thalf_s'])
    return A, indice, fugas


def propagar(A, n0, tiempos_h):
    """N(t) = exp(A·t)·N₀ sobre una malla ORDENADA. Devuelve (n_t, n_especies).

    Se factoriza una exponencial por paso distinto: con malla uniforme es una sola.
    """
    salida = np.zeros((len(tiempos_h), len(n0)))
    cache = {}
    v = np.asarray(n0, dtype=float).copy()
    t_anterior = 0.0
    for k, t in enumerate(tiempos_h):
        dt = t - t_anterior
        clave = round(dt, 9)
        if clave not in cache:
            cache[clave] = expm(A * (dt * 3600.0))
        v = cache[clave] @ v
        salida[k] = v
        t_anterior = t
    return salida


# ---------------------------------------------------------------------------
# 4. Magnitudes derivadas
# ---------------------------------------------------------------------------

def _pico_refinado(t, y):
    """Instante y valor del máximo, con parábola sobre tres nodos."""
    k = int(np.argmax(y))
    if k in (0, len(y) - 1):
        return t[k], y[k]
    y0, y1, y2 = y[k - 1], y[k], y[k + 1]
    denominador = y0 - 2 * y1 + y2
    dk = 0.5 * (y0 - y2) / denominador if denominador else 0.0
    dt = t[1] - t[0]
    return t[k] + dk * dt, y1 - 0.25 * (y0 - y2) * dk


def _cruce(t, p, umbral=UMBRAL_PUREZA):
    """Primer instante en que P supera el umbral, por interpolación lineal."""
    if p[0] >= umbral:
        return 0.0
    i = int(np.argmax(p >= umbral))
    if p[i] < umbral:
        return math.nan
    dt = t[1] - t[0]
    return t[i - 1] + dt * (umbral - p[i - 1]) / (p[i] - p[i - 1])


# ---------------------------------------------------------------------------
# 5. Programa
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description='Propagación analítica del enfriamiento desde los nodos de '
                    'irradiación de una simulación de referencia de ACAB.')
    ap.add_argument('referencia', help='directorio con inp.5, fort.6 y DECAY.dat')
    ap.add_argument('--csv', help='vuelca la tabla del cálculo paramétrico a un CSV')
    ap.add_argument('--dt', type=float, default=0.05,
                    help='paso de la malla fina de enfriamiento, en horas (0,05)')
    ap.add_argument('--tmax', type=float, default=400.0,
                    help='horizonte de la malla fina, en horas (400)')
    args = ap.parse_args()

    ruta_inp5 = os.path.join(args.referencia, 'inp.5')
    ruta_fort6 = os.path.join(args.referencia, 'fort.6')
    ruta_decay = os.path.join(args.referencia, 'DECAY.dat')
    for r in (ruta_inp5, ruta_fort6, ruta_decay):
        if not os.path.exists(r):
            sys.exit('falta %s' % r)

    libreria = leer_decay(ruta_decay)
    especies = especies_TeIXe(libreria)
    nombres = [nombre(n) for n in especies]
    A, indice, fugas = matriz_decaimiento(libreria, especies)
    lam = np.array([libreria[n]['lam'] for n in especies])

    idx_131 = nombres.index('I131')
    yodos = [i for i, s in enumerate(nombres) if s.startswith('I1')]
    masa = np.array([MASA_ATOMICA_YODO.get(s, float(libreria[especies[i]]['A']))
                     for i, s in enumerate(nombres)])
    techo = lam[idx_131] * N_AVOGADRO / MASA_ATOMICA_YODO['I131'] / 1e6   # MBq/g

    print('=' * 74)
    print('SISTEMA')
    print('=' * 74)
    print('  especies Te-I-Xe, A = 127 a 136, con isómeros : %d' % len(especies))
    print('  ramas que salen del conjunto (todas Xe -> Cs) : %d' % len(fugas))
    for k, (f, th) in sorted(fugas.items()):
        print('      %-7s rama %.5f   T1/2 = %.4e s' % (k, f, th))
    print('  techo sin portador (masa atómica, F19)        : %.6e MBq/g' % techo)
    print('  techo con número másico, NO usado             : %.6e MBq/g'
          % (lam[idx_131] * N_AVOGADRO / 131 / 1e6))

    # ---------------- inventarios de irradiación -------------------------
    t_irr, inventario = tabla_atomos(ruta_fort6, ruta_inp5, set(nombres))
    t_irr = np.array(t_irr)
    N0 = np.array([[inventario[k][j] if k in inventario else 0.0
                    for j in range(len(t_irr))] for k in nombres])

    t_enf, actividad_acab = tabla_actividad(ruta_fort6, ruta_inp5, set(nombres))
    t_enf = np.array(t_enf)
    print('  nodos de irradiación (0 incluido)             : %d, hasta %.1f h'
          % (len(t_irr), t_irr[-1]))
    print('  nodos de enfriamiento de ACAB                 : %d, hasta %.1f h'
          % (len(t_enf) - 1, t_enf[-1]))

    # ---------------- 1. verificación contra ACAB ------------------------
    print()
    print('=' * 74)
    print('1. VERIFICACIÓN contra los nodos de enfriamiento que ACAB sí imprime')
    print('=' * 74)
    n_eoi = N0[:, -1]
    Y = np.vstack([n_eoi, propagar(A, n_eoi, t_enf[1:])])
    ACT = Y * lam

    total = positivos = 0
    peor, peor_especie = 0.0, None
    for k in nombres:
        if k not in actividad_acab:
            continue
        a_acab = np.array(actividad_acab[k])
        a_calc = ACT[:, nombres.index(k)]
        m = (a_acab > 1.0) & (a_calc > 1.0)   # por debajo de 1 des/s no hay cifras útiles
        if not m.any():
            continue
        d = 100 * (a_calc[m] / a_acab[m] - 1)
        total += int(m.sum())
        positivos += int((d > 1e-9).sum())
        if abs(d).max() > abs(peor):
            peor, peor_especie = d[np.argmax(abs(d))], k
    print('  comparaciones (especies con A > 1 des/s x nodos) : %d' % total)
    print('  peor residuo                                     : %+.4f %% (%s)'
          % (peor, peor_especie))
    print('  residuos positivos                               : %d' % positivos)

    a131_acab = np.array(actividad_acab['I131'])
    d131 = 100 * (ACT[:, idx_131] / a131_acab - 1)
    print('  I-131: media %+.4f %% | peor %+.4f %% | fuera de ±0,1 %%: %d de %d'
          % (d131.mean(), d131[np.argmax(abs(d131))],
             int((abs(d131) > 0.1).sum()), len(d131)))
    pendiente = np.polyfit(t_enf, d131, 1)[0] * 100
    print('  I-131: deriva del residuo                        : %+.4f %% / 100 h'
          % pendiente)

    # ---------------- 2. control de rama ---------------------------------
    print()
    print('=' * 74)
    print('2. CONTROL DE RAMA — firma de una rama ausente')
    print('=' * 74)
    A2 = A.copy()
    A2[idx_131, nombres.index('TE131M')] = 0.0     # se suprime la β⁻ del ¹³¹ᵐTe
    Y2 = np.vstack([n_eoi, propagar(A2, n_eoi, t_enf[1:])])
    deriva = 100 * (Y2[:, idx_131] / Y[:, idx_131] - 1)
    print('  suprimiendo la β⁻ del Te-131m (79 % de sus desintegraciones):')
    for h in (24, 72, 144):
        k = int(np.argmin(abs(t_enf - h)))
        print('      a %5.1f h : %+.3f %%' % (t_enf[k], deriva[k]))
    print('  una rama ausente NO da un residuo mayor: da un residuo CRECIENTE.')
    k24 = int(np.argmin(abs(t_enf - 24)))
    print('  pendiente entre 24 y 144 h: %+.2f %% / 100 h  (el estimador del capítulo;'
          % (100 * (deriva[-1] - deriva[k24]) / (t_enf[-1] - t_enf[k24])))
    print('   la curva satura, así que un ajuste sobre los 144 h da %+.2f %% / 100 h)'
          % (np.polyfit(t_enf, deriva, 1)[0] * 100))

    # ---------------- 3. el cálculo paramétrico --------------------------------------
    print()
    print('=' * 74)
    print('3. EL CÁLCULO PARAMÉTRICO — los 79 nodos de irradiación, propagados')
    print('=' * 74)
    nt = int(round(args.tmax / args.dt)) + 1
    t = np.arange(nt) * args.dt
    E = expm(A * (args.dt * 3600.0))

    INV = np.zeros((nt, len(especies), len(t_irr)))
    v = N0.copy()
    INV[0] = v
    for k in range(1, nt):
        v = E @ v
        INV[k] = v
    ACTF = INV * lam[None, :, None]

    A131 = ACTF[:, idx_131, :]
    A_yodo = ACTF[:, yodos, :].sum(axis=1)
    P = np.divide(A131, A_yodo, out=np.zeros_like(A131), where=A_yodo > 0)
    n_yodo = INV[:, yodos, :].sum(axis=1)
    frac_iso = np.divide(INV[:, idx_131, :], n_yodo,
                         out=np.zeros_like(A131), where=n_yodo > 0)
    masa_yodo_g = (INV[:, yodos, :] * masa[yodos][None, :, None]).sum(axis=1) / N_AVOGADRO
    A_esp = np.divide(A131, masa_yodo_g, out=np.zeros_like(A131), where=masa_yodo_g > 0)

    filas = []
    for j in range(len(t_irr)):
        if A131[:, j].max() <= 0:          # el nodo t = 0 no es un punto del cálculo paramétrico
            continue
        tm, am = _pico_refinado(t, A131[:, j])
        tc = _cruce(t, P[:, j])
        topt = max(tm, tc)
        aopt = np.interp(topt, t, A131[:, j])
        filas.append((t_irr[j], tm, am, tc, tc - tm, topt, aopt / am,
                      1 - P[0, j], frac_iso[0, j],
                      np.interp(topt, t, frac_iso[:, j]),
                      A_esp[0, j] / 1e6))
    T = np.array(filas)

    cab = ('t_irr', 't_max', 't_cruce', 'separ.', 'A_opt/A_max', '1-P(EOI)',
           'iso(EOI)', 'iso(opt)', 'A_esp(EOI)')
    print('  %7s %8s %9s %9s %12s %11s %9s %9s %12s' % cab)
    for f in T:
        if f[0] in (0.25, 1, 2.5, 5, 10, 20, 35, 50, 65, 80, 95, 110, 125, 140, 150):
            print('  %7.2f %8.3f %9.3f %9.3f %12.4f %11.4e %9.4f %9.4f %12.4e'
                  % (f[0], f[1], f[3], f[4], f[6], f[7], f[8], f[9], f[10]))

    print()
    print('  MONOTONÍA, sobre los %d nodos:' % len(T))
    no_mon_P = sum(1 for j in range(len(t_irr)) if A131[:, j].max() > 0
                   and np.any(np.diff(P[:, j]) < -1e-12))
    max_int = sum(1 for j in range(len(t_irr)) if A131[:, j].max() > 0
                  and int(np.argmax(frac_iso[:, j])) != 0)
    print('      nodos con P(t) no creciente en algún tramo : %d' % no_mon_P)
    print('      nodos con máximo interior de pureza isotópica : %d' % max_int)
    print('      -> ninguna espera recupera la pureza isotópica, para ningún t_irr')
    print()
    print('  t_max recorre %.3f -> %.3f h  (rango %.2f h)'
          % (T[0, 1], T[-1, 1], T[0, 1] - T[-1, 1]))
    print('  t_cruce recorre %.3f -> %.3f h  (rango %.2f h)'
          % (T[0, 3], T[-1, 3], T[-1, 3] - T[0, 3]))
    rec_cruce = T[-1, 3] - T[0, 3]
    rec_max = T[0, 1] - T[-1, 1]
    print('  separación entre t_max y t_cruce a %.0f h: %.2f h' % (T[-1, 0], T[-1, 4]))
    print('  reparto de esa separación: el t_cruce pone el %.0f %% (%.2f h de recorrido)'
          % (100 * rec_cruce / (rec_cruce + rec_max), rec_cruce))
    print('                             y el adelanto del máximo el %.0f %% (%.2f h)'
          % (100 * rec_max / (rec_cruce + rec_max), rec_max))
    print('  precio de la espera a %.0f h de irradiación: %.2f %% de la actividad de pico'
          % (T[-1, 0], 100 * (1 - T[-1, 6])))

    # ---------------- 4. control asintótico (E3-bis) ---------------------
    print()
    print('=' * 74)
    print('4. CONTROL ASINTÓTICO (E3-bis) — verificación, no resultado')
    print('=' * 74)
    print('  Congelada la masa de yodo, A_esp debe decaer exactamente como el I-131.')
    DT_L, NT_L = 10.0, 6001
    tl = np.arange(NT_L) * DT_L
    EL = expm(A * (DT_L * 3600.0))
    v = n_eoi.copy()
    a131_l = np.zeros(NT_L)
    masa_l = np.zeros(NT_L)
    for k in range(NT_L):
        if k:
            v = EL @ v
        a131_l[k] = v[idx_131] * lam[idx_131]
        masa_l[k] = (v[yodos] * masa[yodos]).sum()
    aesp_l = a131_l / masa_l
    lam_h = lam[idx_131] * 3600.0
    for t0, t1 in ((1200, 2400), (2400, 4800), (4800, 8000)):
        m = (tl >= t0) & (tl <= t1)
        p = np.polyfit(tl[m], np.log(aesp_l[m]), 1)
        print('      %5d-%5d h : cociente con lambda(I-131) = %.5f' % (t0, t1, -p[0] / lam_h))
    kmin = int(np.argmin(masa_l))
    print('  la masa de yodo NO se congela en seguida: mínimo en %.0f h al %.5f'
          % (tl[kmin], masa_l[kmin] / masa_l[0]))
    print('  y se estabiliza en %.5f hacia las %.0f h (manda el Te-127m, 106,1 d)'
          % (masa_l[-1] / masa_l[0], tl[-1]))
    print('  medir este control a horizonte corto lo da por fallado sin que falle nada.')

    if args.csv:
        cabecera = ('t_irr_h,t_max_h,A_max_des_s,t_cruce_h,separacion_h,t_opt_h,'
                    'A_opt_sobre_A_max,uno_menos_P_EOI,frac_iso_EOI,frac_iso_opt,'
                    'A_esp_EOI_MBq_g')
        np.savetxt(args.csv, T, delimiter=',', header=cabecera, comments='',
                   fmt='%.8g')
        print()
        print('  tabla de los %d nodos volcada en %s' % (len(T), args.csv))


if __name__ == '__main__':
    main()
