#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
probe_collaps_mesh.py
=====================
Sondea la malla de energia que COLLAPS usa INTERNAMENTE para ILIB=IESF=12,
usando la opcion ISTOP=1 (manual ACAB v.2008, COLL.inp tarjeta 9: "the code only
writes the file FLUX.inf"). Cada sondeo es una ejecucion de segundos que no
colapsa nada: solo imprime el diagnostico del espectro.

IDEA
----
Con flujo 1,0 en un unico grupo g y 0 en el resto, la energia media que declara
FLUX.inf es la energia representativa que COLLAPS asigna a ese grupo, dividida
por el denominador que use (FLUX.inf imprime tambien FTL, asi que el denominador
se despeja). Repetido grupo a grupo, se obtiene la malla LEIDA DEL CODIGO, que es
la unica fuente que no depende ni de la reconstruccion 1/v ni de una tabla externa.

MODOS
-----
  gen   : escribe los COLL.inp de sondeo en un directorio
  parse : lee los FLUX.inf resultantes y los contrasta con la malla de referencia

USO
---
  # 1) generar (por defecto, 12 sondeos representativos; --all para los 211)
  python3 probe_collaps_mesh.py gen  probes/  [--all]

  # 2) ejecutar cada sondeo (fuera de este script; ejemplo en bash)
  #    for d in probes/p*/ ; do (cd $d && cp ../../XSBL.dat . && collaps) ; done
  #    NOTA: aunque ISTOP=1 no colapsa, COLLAPS abre XSBL.dat para leer la
  #    cabecera (IHEAD lineas); si su instalacion no lo exige, omita la copia.

  # 3) analizar
  python3 probe_collaps_mesh.py parse probes/ vitaminj_plus_211_eV.txt

QUE ESPERAR
-----------
  - Si E_declarada(g) coincide con (Elo+Ehi)/2 de la malla de referencia en los
    211 grupos, la malla queda confirmada CONTRA EL PROPIO CODIGO y el punto 1.2
    se cierra sin depender de ninguna fuente externa.
  - Si coincide con sqrt(Elo*Ehi), lo mismo, con convenio geometrico.
  - Si no coincide con ninguna, la discrepancia esta en la malla y hay que
    investigarla ANTES de dar por bueno el XSECTION.dat.
  - Los sondeos de escala (x2) y de par (dos grupos) sirven para fijar que es
    FTL: si FTL escala x2 con el flujo, es proporcional al total; si no, depende
    de la forma del espectro.
"""

import os
import re
import sys
import numpy as np

NG = 211
ILIB, IESF, IHEAD, FF = 12, 12, 16, 0

# grupos representativos: termicos, epitermicos, resonantes, rapidos y extremos
DEFAULT_GROUPS = [211, 210, 209, 205, 200, 190, 180, 150, 100, 50, 36, 1]


def write_coll_inp(path, phi):
    with open(path, "w") as f:
        f.write("%4d%4d\n" % (ILIB, IESF))
        f.write("%d\n" % IHEAD)
        f.write("0 0 0 0\n")
        f.write("%4d %3d\n" % (-NG, FF))
        for i in range(0, NG, 6):
            f.write("".join("%12.5E" % x for x in phi[i:i + 6]) + "\n")
        f.write("0\n")   # IUNC3G = 0
        f.write("1\n")   # ISTOP  = 1  <-- solo escribe FLUX.inf


def gen(outdir, groups):
    os.makedirs(outdir, exist_ok=True)
    cases = []
    for g in groups:
        phi = np.zeros(NG); phi[g - 1] = 1.0
        cases.append(("p_g%03d" % g, phi))
    # control de escala: mismo grupo, flujo x2 -> FTL proporcional al total?
    phi = np.zeros(NG); phi[179] = 2.0
    cases.append(("p_g180_x2", phi))
    # control de aditividad: dos grupos muy separados
    phi = np.zeros(NG); phi[210] = 1.0; phi[179] = 1.0
    cases.append(("p_g211_g180", phi))
    for name, phi in cases:
        d = os.path.join(outdir, name)
        os.makedirs(d, exist_ok=True)
        write_coll_inp(os.path.join(d, "COLL.inp"), phi)
    print("Escritos %d sondeos en %s" % (len(cases), outdir))
    print("Cada uno con ISTOP=1: COLLAPS solo escribira FLUX.inf.")


def read_flux_inf(path):
    """Extrae (flujo total, energia media en MeV, FTL) de un FLUX.inf."""
    t = open(path, encoding='latin-1').read()
    num = r'([-+]?\d*\.?\d+[EeDd][-+]?\d+)'
    m1 = re.search(r'REAL TOTAL FLUX AND AVERAGE ENERGY.*?' + num + r'\s+' + num, t, re.S)
    m2 = re.search(r'\*\s*FTL\s*\*\s*' + num, t, re.S)
    f = lambda s: float(s.replace('D', 'E').replace('d', 'e'))
    return (f(m1.group(1)), f(m1.group(2)), f(m2.group(1)) if m2 else float('nan'))


def read_mesh(path):
    Ehi, Elo = [], []
    for line in open(path):
        if line.lstrip().startswith('#') or not line.strip():
            continue
        g, hi, lo = line.split()
        Ehi.append(float(hi)); Elo.append(float(lo))
    return np.array(Ehi), np.array(Elo)


def parse(outdir, meshfile):
    Ehi, Elo = read_mesh(meshfile)
    conv = {'aritmetica (Elo+Ehi)/2': 0.5 * (Ehi + Elo),
            'geometrica sqrt(Elo*Ehi)': np.sqrt(Ehi * Elo),
            'letargia (Ehi-Elo)/ln(Ehi/Elo)': (Ehi - Elo) / np.log(Ehi / Elo),
            'frontera superior Ehi': Ehi,
            'frontera inferior Elo': Elo}
    rows, extra = [], []
    for name in sorted(os.listdir(outdir)):
        fi = os.path.join(outdir, name, "FLUX.inf")
        if not os.path.exists(fi):
            continue
        tot, emed, ftl = read_flux_inf(fi)
        m = re.fullmatch(r'p_g(\d{3})', name)
        if m:
            rows.append((int(m.group(1)), tot, emed * 1e6, ftl))   # eV
        else:
            extra.append((name, tot, emed * 1e6, ftl))
    if not rows:
        print("No se encontro ningun FLUX.inf en %s" % outdir); return

    print("=== Energia representativa declarada por COLLAPS, grupo a grupo ===")
    print("%5s %12s %12s %12s %10s" % ("grupo", "E_decl(eV)", "FTL", "flujo tot", "E*FTL(eV)"))
    for g, tot, e, ftl in rows:
        print("%5d %12.5E %12.5E %12.5E %12.5E" % (g, e, ftl, tot, e * ftl / tot))
    print("\n=== Contraste con la malla de referencia (%s) ===" % os.path.basename(meshfile))
    gs = np.array([r[0] for r in rows]) - 1
    for corr, lab in ((1.0, "E_decl"), (None, "E_decl*FTL/flujo")):
        vals = np.array([r[2] if corr else r[2] * r[3] / r[1] for r in rows])
        print("  -- interpretando la salida como %s --" % lab)
        for cname, cvals in conv.items():
            r = vals / cvals[gs]
            print("     %-32s razon media %.4f   dispersion %.2f %%"
                  % (cname, r.mean(), 100 * (r.max() / r.min() - 1)))
    if extra:
        print("\n=== Controles de escala y aditividad ===")
        for name, tot, e, ftl in extra:
            print("  %-14s flujo %.5E  E %.5E eV  FTL %.5E" % (name, tot, e, ftl))
        print("  (p_g180_x2 con FTL doble que p_g180 => FTL proporcional al flujo total)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    if sys.argv[1] == "gen":
        groups = list(range(1, NG + 1)) if "--all" in sys.argv else DEFAULT_GROUPS
        gen(sys.argv[2], groups)
    elif sys.argv[1] == "parse":
        parse(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "vitaminj_plus_211_eV.txt")
    else:
        print(__doc__)
