# run_all.ps1 — lanza ACAB en cada subcarpeta del barrido.
# Ajusta $ACAB_EXE a la ruta del ejecutable de ACAB.
$ACAB_EXE = "acab.exe"
$dirs = @(
  'TeO2x0.049810719266786214',
  'TeO2x0.664474995018928',
  'TeO2x0.8477784419207013',
  'TeO2x1',
  'TeO2x996.2143853357243'
)
foreach ($d in $dirs) {
  Push-Location (Join-Path $PSScriptRoot $d)
  & $ACAB_EXE *> run.log
  Pop-Location
}
