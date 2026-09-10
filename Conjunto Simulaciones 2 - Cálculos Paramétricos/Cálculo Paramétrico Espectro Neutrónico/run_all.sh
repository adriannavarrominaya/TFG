#!/usr/bin/env bash
# run_all.sh — lanza ACAB en cada subcarpeta del barrido.
# Ajusta ACAB_EXE a la ruta del ejecutable de ACAB.
set -u
ACAB_EXE="${ACAB_EXE:-acab}"
cd "$(dirname "$0")"
for d in "TeO269_LANL-OWR" "TeO270_Cf252" "TeO2112_MURR-G1" "TeO2171_HFR-C3" "TeO2172_Phenix" "TeO2175_ITER-DT" "TeO2238_HFIR-VXF3-AD" "TeO2621_SCK-BR2" "TeO2640_LR-0-Void"; do
  (cd "$d" && "$ACAB_EXE" > run.log 2>&1)
done
