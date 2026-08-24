#!/usr/bin/env bash
# run_all.sh — lanza ACAB en cada subcarpeta del barrido.
# Ajusta ACAB_EXE a la ruta del ejecutable de ACAB.
set -u
ACAB_EXE="${ACAB_EXE:-acab}"
cd "$(dirname "$0")"
for d in "TeO2Tirr_0p5h" "TeO2Tirr_0p25h" "TeO2Tirr_1h" "TeO2Tirr_2h" "TeO2Tirr_4h" "TeO2Tirr_7h" "TeO2Tirr_12h" "TeO2Tirr_20h" "TeO2Tirr_32h" "TeO2Tirr_50h" "TeO2Tirr_75h" "TeO2Tirr_100h" "TeO2Tirr_125h" "TeO2Tirr_150h"; do
  (cd "$d" && "$ACAB_EXE" > run.log 2>&1)
done
