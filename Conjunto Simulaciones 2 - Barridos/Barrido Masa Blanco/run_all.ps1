# run_all.ps1 — lanza ACAB en cada subcarpeta del barrido.
# Ajusta $ACAB_EXE a la ruta del ejecutable de ACAB.
$ACAB_EXE = "acab.exe"
$dirs = @(
  'TeO2m0.123g',
  'TeO2m1.000g',
  'TeO2m10.000g',
  'TeO2m100.000g',
  'TeO2m500.000g',
  'TeO2m1000.000g'
)
foreach ($d in $dirs) {
  Push-Location (Join-Path $PSScriptRoot $d)
  & $ACAB_EXE *> run.log
  Pop-Location
}
