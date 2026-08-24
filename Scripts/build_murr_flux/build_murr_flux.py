#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_murr_flux.py  --  version 2 (malla real)
==============================================
Construye el vector de flujo neutronico del MURR (Exp. 1 de Haffner et al.) en la
malla VITAMIN-J+ de 211 grupos y genera el COLL.inp listo para COLLAPS.

CAMBIO FRENTE A LA VERSION 1
----------------------------
La v.1 RECONSTRUIA las fronteras de energia invirtiendo el comportamiento 1/v del
H-1. La v.2 LEE la malla real de un fichero con procedencia citable
(`vitaminj_plus_211_eV.txt`, de ENEA ADPFISS-LP1-107 Tabla 2, validada contra la
Tabla III.1 del manual ACAB v.2008). El comportamiento 1/v del H-1 se conserva,
pero degradado a CONTROL de la malla, no a fuente de la malla.

Efecto medido del cambio sobre el espectro del Exp. 1: < 0,11 % en todos los
canales probados (1/v y resonantes). La v.1 no sesgaba este caso; simplemente no
era citable, y era ciega por encima de 0,297 MeV.

FUENTES DE CADA VALOR
---------------------
  Malla        : fichero externo (ver su cabecera).
  phi_termico  : 6,5e13 n/cm2/s   -- Haffner, Miller & Morris,
  phi_epiterm. : 1,7e12 n/cm2/s      Appl. Radiat. Isot. 151 (2019) 52-61, Exp. 1.
  kT           : 0,0253 eV (energia termica de referencia a 2200 m/s). Con el
                 kT de la temperatura de la libreria (293,16 K -> 0,025263 eV)
                 sigma(Te-130) cambia +0,03 % : irrelevante. Declarado, no ajustado.
  EPI_LO       : 0,5 eV -- corte de cadmio, convenio del reparto termico/epitermico.
  EPI_HI       : 1,0e5 eV \\ HIPOTESIS DEL MODELO, no dato del paper. Variar EPI_HI
  rapido = 0   :          / entre 5e4 y 2e5 eV mueve sigma(Te-130) < 0,02 %.
                 Consecuencia declarable: este espectro excluye las reacciones umbral.

USO
---
   python3 build_murr_flux.py XSBL.dat [vitaminj_plus_211_eV.txt]
Genera:
   flux_murr_211.txt  -- vector de flujo, 6 por linea, orden descendente
   COLL.inp           -- fichero de entrada completo de COLLAPS
"""

import sys
import hashlib
import numpy as np

NG       = 211
kT       = 0.0253      # eV
PHI_TH   = 6.5e13      # n/cm2/s
PHI_EPI  = 1.7e12      # n/cm2/s
EPI_LO   = 0.5         # eV
EPI_HI   = 1.0e5       # eV

# Tarjetas de COLL.inp (manual ACAB v.2008, cap. III)
ILIB, IESF = 12, 12    # Vitamin-J+ (211 grupos) para libreria y para flujo
IHEAD      = 16        # lineas de cabecera de XSBL.dat
FF         = 0         # flujo escalar total por grupo [n/cm2-s]

KEY_H1   = '10010 1020'    # H-1  (n,g) H-2      -> control 1/v de la malla
KEY_TE_G = '521300 1020'   # Te-130 (n,g) Te-131
KEY_TE_M = '521300 1021'   # Te-130 (n,g) Te-131m


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def read_mesh(path):
    """Lee la malla real: devuelve Ehi, Elo en eV, grupos 1..211 (descendente)."""
    Ehi, Elo = [], []
    for line in open(path):
        if line.lstrip().startswith('#') or not line.strip():
            continue
        g, hi, lo = line.split()
        assert int(g) == len(Ehi) + 1, "malla desordenada en el grupo %s" % g
        Ehi.append(float(hi)); Elo.append(float(lo))
    assert len(Ehi) == NG, "la malla tiene %d grupos, se esperaban %d" % (len(Ehi), NG)
    Ehi, Elo = np.array(Ehi), np.array(Elo)
    assert np.all(Ehi > Elo) and np.all(np.diff(Ehi) < 0), "malla no monotona"
    assert np.allclose(Ehi[1:], Elo[:-1]), "fronteras no contiguas"
    return Ehi, Elo


def read_multigroup(path, key, ng=NG):
    """Lee los 'ng' valores multigrupo de una reaccion de XSBL.dat (formato EAF)."""
    lines = open(path, encoding='latin-1').read().splitlines()
    k = key.replace(' ', '')
    for i, l in enumerate(lines):
        if l.replace(' ', '').startswith(k):
            vals, j = [], i + 1
            while len(vals) < ng and j < len(lines):
                toks = lines[j].split()
                ok, nums = len(toks) > 0, []
                for t in toks:
                    try:
                        nums.append(float(t))
                    except ValueError:
                        ok = False
                        break
                if ok:
                    vals.extend(nums)
                elif vals:
                    break
                j += 1
            return np.array(vals[:ng])
    raise KeyError("No se encontro la reaccion %r en %s" % (key, path))


def check_mesh_against_library(H, Ehi, Elo):
    """CONTROL (no fuente): el H-1(n,g) es un 1/v puro sin resonancias, luego su
    valor medio con peso plano sobre [Elo,Ehi] vale
         sigma_avg = 2*sigma0*sqrt(E0)/(sqrt(Ehi)+sqrt(Elo)).
    Se ajusta sigma0 a la libreria en la region 1/v y se comprueba (a) que el
    ajuste es consistente grupo a grupo y (b) que el sigma0 resultante coincide
    con el valor aceptado del H-1 a 2200 m/s (0,3326 b). Ademas se comprueba que
    la frontera de 20 MeV de la malla cae donde la libreria deja de tener datos."""
    E0 = 0.0253
    sel = np.arange(184, NG)               # grupos 185..211: region 1/v limpia
    s0 = np.mean(H[sel] * (np.sqrt(Ehi[sel]) + np.sqrt(Elo[sel])) / (2 * np.sqrt(E0)))
    pred = 2 * s0 * np.sqrt(E0) / (np.sqrt(Ehi) + np.sqrt(Elo))
    dev = np.max(np.abs(H[sel] / pred[sel] - 1)) * 100
    # frontera de 20 MeV: la libreria es nula por encima
    first_nz = int(np.argmax(H > 0)) + 1
    print("=== CONTROL 1: malla real contra el H-1(n,g) de la libreria ===")
    print("  sigma0 ajustado en los grupos 185-211 : %.4f b  (aceptado 0,3326 b, dif %+.2f %%)"
          % (s0, 100 * (s0 / 0.3326 - 1)))
    print("  dispersion grupo a grupo del ajuste   : %.3f %%" % dev)
    print("  primer grupo con datos en la libreria : %d  -> E_sup = %.4g MeV (frontera de 20 MeV)"
          % (first_nz, Ehi[first_nz - 1] * 1e-6))
    return s0, dev


def build_flux(Ehi, Elo):
    """Flujo integrado por grupo [n/cm2-s]: maxwelliana de flujo + 1/E + rapido nulo."""
    # Maxwelliana de flujo phi(E) ~ E*exp(-E/kT);  INT = -kT*(E+kT)*exp(-E/kT)
    F = lambda E: -kT * (E + kT) * np.exp(-E / kT)
    th_w = F(Ehi) - F(Elo)
    # 1/E recortada a [EPI_LO, EPI_HI]: peso = ln(Ehi/Elo) del solape
    lo = np.maximum(Elo, EPI_LO)
    hi = np.minimum(Ehi, EPI_HI)
    ep_w = np.where(hi > lo, np.log(np.where(hi > lo, hi / np.maximum(lo, 1e-30), 1.0)), 0.0)
    return PHI_TH * th_w / th_w.sum() + PHI_EPI * ep_w / ep_w.sum()


def write_coll_inp(phi, path="COLL.inp"):
    """COLL.inp completo. Tarjetas segun el manual ACAB v.2008, cap. III."""
    with open(path, "w") as f:
        f.write("%4d%4d\n" % (ILIB, IESF))                  # Card #1
        f.write("%d\n" % IHEAD)                             # Card #2
        f.write("0 0 0 0\n")                                # Card #3 ISFIS=0
        f.write("%4d %3d\n" % (-NG, FF))                    # Card #5 NGROUP<0, FF=0
        for i in range(0, NG, 6):                           # Card #7 (6E12.5)
            f.write("".join("%12.5E" % x for x in phi[i:i + 6]) + "\n")
        f.write("0\n")                                      # Card #8 IUNC3G=0
        f.write("0\n")                                      # Card #9 ISTOP=0


def main():
    xsbl = sys.argv[1] if len(sys.argv) > 1 else "XSBL.dat"
    mesh = sys.argv[2] if len(sys.argv) > 2 else "vitaminj_plus_211_eV.txt"

    print("Libreria : %s\n           sha256 %s" % (xsbl, sha256(xsbl)))
    print("Malla    : %s\n           sha256 %s\n" % (mesh, sha256(mesh)))

    Ehi, Elo = read_mesh(mesh)
    H    = read_multigroup(xsbl, KEY_H1)
    Te_g = read_multigroup(xsbl, KEY_TE_G)
    Te_m = read_multigroup(xsbl, KEY_TE_M)

    check_mesh_against_library(H, Ehi, Elo)

    phi = build_flux(Ehi, Elo)
    den = phi.sum()
    Emed = (phi * 0.5 * (Ehi + Elo)).sum() / den

    print("\n=== Espectro MURR construido (Exp. 1) ===")
    print("  flujo total               = %.4E n/cm2/s   (paper: 6,5E13 + 1,7E12 = 6,67E13)" % den)
    print("  energia media (E_centro)  = %.4E eV" % Emed)
    print("  fraccion en el grupo 211  = %.3f %%   [1E-5 ; 0,1 eV]" % (100 * phi[210] / den))
    print("  fraccion en 205-211       = %.3f %%   [1E-5 ; 1,445 eV]" % (100 * phi[204:].sum() / den))
    print("  fraccion en 190-211       = %.3f %%   [1E-5 ; 61,44 eV]" % (100 * phi[189:].sum() / den))

    sig_g = (Te_g * phi).sum() / den
    sig_m = (Te_m * phi).sum() / den
    print("\n=== CONTROL 2: colapso a 1 grupo del Te-130(n,gamma) ===")
    print("  ->Te-131   = %.6f b" % sig_g)
    print("  ->Te-131m  = %.7f b" % sig_m)
    print("  TOTAL      = %.6f b" % (sig_g + sig_m))
    print("  equivale a sigma_0 = TOTAL/(sqrt(pi)/2) = %.4f b" % ((sig_g + sig_m) / 0.8862269))
    print("  (Tabla 3 de Haffner et al.: NNDC 2014 0,195 ; JEFF-3.2 0,29 ; Firestone 0,29)")

    with open("flux_murr_211.txt", "w") as fo:
        fo.write("# Flujo MURR Exp.1, malla VITAMIN-J+ 211 grupos, orden descendente.\n")
        fo.write("# COLLAPS: ILIB=IESF=12, NGROUP=-211, FF=0 [n/cm2-s por grupo].\n")
        for i in range(0, NG, 6):
            fo.write(" ".join("%.5E" % x for x in phi[i:i + 6]) + "\n")
    write_coll_inp(phi)
    print("\nEscritos: flux_murr_211.txt  y  COLL.inp")


if __name__ == "__main__":
    main()
