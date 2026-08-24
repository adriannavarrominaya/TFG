#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CÓDIGO DE APOYO — ejemplo comentado
===================================

Esto NO es la suite. Es el tipo de comprobación independiente que se ha usado
para verificar ACAB: unas pocas decenas de líneas que resuelven el mismo
problema por otro camino y contrastan el resultado.

Su valor está precisamente en ser corto y en no compartir nada con la
herramienta: lee los mismos ficheros de entrada que ACAB y resuelve la cadena
en forma cerrada.

Tres ejemplos, de menor a mayor:

  1. Leer los datos nucleares de los ficheros de entrada.
  2. Resolver la cadena A = 131 analíticamente y contrastarla con ACAB.
  3. Comprobar el balance de la cadena: lo que hay debe ser lo producido.

Uso:  python apoyo_ejemplo.py <carpeta_de_la_simulacion>
"""

import sys
import math
import numpy as np
from scipy.linalg import expm

LN2 = math.log(2.0)
BARN = 1.0e-24            # cm2
NA = 6.02214076e23        # 1/mol


# ---------------------------------------------------------------------------
# 1. LEER LOS DATOS NUCLEARES DE LOS FICHEROS DE ENTRADA
# ---------------------------------------------------------------------------
# La idea es no fiarse de ningún valor tabulado a mano: las secciones eficaces
# y las semividas se leen de los mismos ficheros que ACAB va a usar, de modo
# que si la biblioteca cambia, la comprobación cambia con ella.

def seccion(carpeta, zaid, mt):
    """Sección eficaz colapsada, en barnios, del `XSECTION.dat`.

    zaid: identificador del nucleido, p. ej. 521300 para el Te-130.
    mt:   canal, 1020 = (n,gamma) al fundamental, 1021 = al metaestable.
    """
    lineas = open(f'{carpeta}/XSECTION.dat', 'rb').read().decode('latin-1')
    lineas = lineas.replace('\r', '').split('\n')
    for i, l in enumerate(lineas):
        if l[:8].strip() == str(zaid) and l[8:13].strip() == str(mt):
            # el valor está en una de las líneas siguientes
            for k in range(i + 1, i + 6):
                try:
                    return float(lineas[k].strip())
                except ValueError:
                    continue
    return None                      # el canal no existe en esta biblioteca


def decaimiento(carpeta, zaid):
    """Devuelve (semivida en s, rama de transición isomérica) del `DECAY.dat`."""
    lineas = open(f'{carpeta}/DECAY.dat', 'rb').read().decode('latin-1')
    lineas = lineas.replace('\r', '').split('\n')
    for l in lineas:
        if len(l) > 12 and l[6:12].strip() == str(zaid):
            campos = l[12:].split()
            return float(campos[1]), float(campos[6])
    return None, None


# ---------------------------------------------------------------------------
# 2. RESOLVER LA CADENA ANALÍTICAMENTE Y CONTRASTAR CON ACAB
# ---------------------------------------------------------------------------
# La cadena del I-131 es lineal, así que se puede escribir como  y' = M y + b
# y resolver en forma cerrada:   y(t) = e^{Mt} y0 + M^-1 (e^{Mt} - I) b
#
# Nucleidos:  0 = Te-131,  1 = Te-131m,  2 = I-131
#
# Durante la irradiación hay término fuente (b) y quemado por captura;
# durante el enfriamiento no hay flujo, así que b = 0 y desaparece el quemado.
# Esa distinción es importante: aplicar el quemado también al enfriamiento es
# un error fácil de cometer y da un desacuerdo del orden del 0,2 %.

def cadena_i131(carpeta, flujo, n_te130):
    """Construye las matrices de irradiación y enfriamiento y el término fuente."""
    sigma_g = seccion(carpeta, 521300, 1020)     # Te-130 -> Te-131
    sigma_m = seccion(carpeta, 521300, 1021)     # Te-130 -> Te-131m

    t_g, _ = decaimiento(carpeta, 521310)        # Te-131
    t_m, rama_it = decaimiento(carpeta, 521311)  # Te-131m, con su rama IT
    t_i, _ = decaimiento(carpeta, 531310)        # I-131

    lam_g, lam_m, lam_i = LN2 / t_g, LN2 / t_m, LN2 / t_i

    # Destrucción por captura durante la irradiación. Ojo: algunos canales no
    # existen en todas las bibliotecas, y `seccion` devuelve None en ese caso.
    def sigma_o_cero(zaid):
        s = seccion(carpeta, zaid, 1020)
        return 0.0 if s is None else s * BARN * flujo

    d_g = sigma_o_cero(521310)
    d_m = sigma_o_cero(521311)
    d_i = sigma_o_cero(531310)

    # Matriz durante la irradiación (con quemado)
    M_irr = np.array([
        [-(lam_g + d_g),  rama_it * lam_m,        0.0],
        [ 0.0,           -(lam_m + d_m),          0.0],
        [ lam_g,         (1 - rama_it) * lam_m,  -(lam_i + d_i)],
    ])

    # Matriz durante el enfriamiento (sin flujo: sin quemado)
    M_enf = np.array([
        [-lam_g,          rama_it * lam_m,        0.0],
        [ 0.0,           -lam_m,                  0.0],
        [ lam_g,         (1 - rama_it) * lam_m,  -lam_i],
    ])

    # Término fuente: las capturas sobre el Te-130, que es prácticamente constante
    fuente = np.array([n_te130 * sigma_g * BARN * flujo,
                       n_te130 * sigma_m * BARN * flujo,
                       0.0])

    return M_irr, M_enf, fuente, lam_i


def actividad_i131(carpeta, flujo, n_te130, t_irr_h, t_enf_h):
    """Actividad del I-131 en Bq/cm3, t_enf_h horas después de apagar."""
    M_irr, M_enf, b, lam_i = cadena_i131(carpeta, flujo, n_te130)

    # Irradiación: solución con término fuente, partiendo de inventario nulo
    E = expm(M_irr * (t_irr_h * 3600.0))
    y_eoi = np.linalg.inv(M_irr) @ (E - np.eye(3)) @ b

    # Enfriamiento: decaimiento puro desde el inventario del apagado
    y = expm(M_enf * (t_enf_h * 3600.0)) @ y_eoi

    return lam_i * y[2]


# ---------------------------------------------------------------------------
# 3. EL BALANCE DE LA CADENA
# ---------------------------------------------------------------------------
# Comprobación independiente de la anterior, y más robusta porque no depende
# de resolver ninguna ecuación: todo átomo de la cadena A = 131 que existe
# tuvo que producirse por una captura sobre el Te-130. Basta sumar el
# inventario y compararlo con el número de capturas.
#
# Es la comprobación que destapó la dependencia de la malla: cuando el
# cociente se aparta del 1, hay átomos que ACAB no ha producido.

def balance_cadena(inventario_eoi, n_te130, sigma_total, flujo, t_irr_h):
    """Cociente entre lo que hay y lo que debería haberse producido.

    inventario_eoi: dict {nucleido: átomos} con TODA la cadena A = 131,
                    incluidos los sumideros estables (Xe-131) y lo que se ha
                    ido a A = 132 por captura. Si se omite alguno, el balance
                    no cierra y parece un error del código cuando no lo es.
    """
    hay = sum(inventario_eoi.values())
    producido = n_te130 * sigma_total * BARN * flujo * (t_irr_h * 3600.0)
    return hay / producido


# ---------------------------------------------------------------------------
# Ejemplo de uso
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    carpeta = sys.argv[1] if len(sys.argv) > 1 else '.'

    # Estos tres datos se leen del `inp.5`: flujo del bloque 3 y densidad
    # atómica del bloque 5, multiplicada por la abundancia del Te-130.
    flujo = 1.0038e14          # n/cm2·s
    n_te = 4.6448e-4 * 1e24    # átomos de Te por cm3
    _, _, abundancia = None, None, 0.33799
    n_te130 = n_te * abundancia

    print(f'Datos leídos de {carpeta}:')
    print(f'  sigma(Te-130 -> Te-131 )  = {seccion(carpeta, 521300, 1020)} b')
    print(f'  sigma(Te-130 -> Te-131m)  = {seccion(carpeta, 521300, 1021)} b')
    print(f'  T1/2(Te-131m)             = {decaimiento(carpeta, 521311)[0]} s')
    print()

    for t_enf in (0.0, 2.0, 24.0):
        a = actividad_i131(carpeta, flujo, n_te130, t_irr_h=150.0, t_enf_h=t_enf)
        print(f'  A(I-131) a EOI + {t_enf:5.1f} h = {a / 1e6 / 0.12317:12.4f} MBq/g')

    print()
    print('Contrastar estos valores con la tabla del `fort.6`.')
    print('Acuerdo esperado: del orden del 0,04 %, que es el suelo de impresión.')
