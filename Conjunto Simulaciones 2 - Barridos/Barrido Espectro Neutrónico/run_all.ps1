# run_all.ps1 — lanza ACAB en cada subcarpeta del barrido.
# Ajusta $ACAB_EXE a la ruta del ejecutable de ACAB.
$ACAB_EXE = "acab.exe"
$dirs = @(
  'TeO269_LANL-OWR',
  'TeO270_Cf252',
  'TeO2112_MURR-G1',
  'TeO2171_HFR-C3',
  'TeO2172_Phenix',
  'TeO2175_ITER-DT',
  'TeO2238_HFIR-VXF3-AD',
  'TeO2621_SCK-BR2',
  'TeO2640_LR-0-Void'
)
foreach ($d in $dirs) {
  Push-Location (Join-Path $PSScriptRoot $d)
  & $ACAB_EXE *> run.log
  Pop-Location
}
