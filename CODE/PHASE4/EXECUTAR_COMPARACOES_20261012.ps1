param([ValidateSet('Preflight','Smoke','Full','Summary')][string]$Mode = 'Preflight')
$ErrorActionPreference = 'Stop'
$pythonExe = 'C:\Users\Admin\anaconda3\envs\spyder_dl8\python.exe'
$scriptPath = Join-Path $PSScriptRoot 'comparacoes_otimizacao_20261012.py'
if (-not (Test-Path -LiteralPath $pythonExe)) { throw "Python ausente: $pythonExe" }
$logDir = Join-Path $PSScriptRoot 'comparacoes_otimizacao_20261012\console_logs'
New-Item -ItemType Directory -Path $logDir -Force | Out-Null
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$logPath = Join-Path $logDir "$($Mode)_$stamp.log"
$env:PYTHONUNBUFFERED = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:MPLBACKEND = 'Agg'
$runArgs = @($scriptPath)
switch ($Mode) {
    'Preflight' { $runArgs += 'preflight' }
    'Smoke' { $runArgs += @('run','--smoke') }
    'Full' {
        $smokeComplete = Join-Path $PSScriptRoot 'comparacoes_otimizacao_20261012\smoke\study_complete.json'
        if (-not (Test-Path -LiteralPath $smokeComplete)) { throw 'Execute -Mode Smoke com sucesso antes de Full.' }
        $runArgs += 'run'
    }
    'Summary' { $runArgs += 'summarize' }
}
Write-Host "Console em tempo real. Log: $logPath"
# PowerShell 5.1 may promote native stderr warnings to errors. Exit code is authoritative.
$ErrorActionPreference = 'Continue'
& $pythonExe -u @runArgs 2>&1 | Tee-Object -FilePath $logPath
$runExit = $LASTEXITCODE
$ErrorActionPreference = 'Stop'
if ($runExit -ne 0) { throw "Execucao falhou (exit=$runExit). Consulte $logPath" }
Write-Host "CONCLUIDO $Mode. Log: $logPath"
