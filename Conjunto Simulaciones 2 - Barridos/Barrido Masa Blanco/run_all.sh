#!/usr/bin/env bash
# run_all.sh — lanza ACAB en cada subcarpeta del barrido.
# Ajusta ACAB_EXE a la ruta del ejecutable de ACAB.
set -u
ACAB_EXE="${ACAB_EXE:-acab}"
cd "$(dirname "$0")"
for d in "TeO2m0.123g" "TeO2m1.000g" "TeO2m10.000g" "TeO2m100.000g" "TeO2m500.000g" "TeO2m1000.000g"; do
  (cd "$d" && "$ACAB_EXE" > run.log 2>&1)
done
