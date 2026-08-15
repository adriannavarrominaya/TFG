#!/usr/bin/env bash
# run_all.sh — lanza ACAB en cada subcarpeta del barrido.
# Ajusta ACAB_EXE a la ruta del ejecutable de ACAB.
set -u
ACAB_EXE="${ACAB_EXE:-acab}"
cd "$(dirname "$0")"
for d in "TeO2x0.049810719266786214" "TeO2x0.664474995018928" "TeO2x0.8477784419207013" "TeO2x1" "TeO2x996.2143853357243"; do
  (cd "$d" && "$ACAB_EXE" > run.log 2>&1)
done
