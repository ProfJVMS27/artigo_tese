[CmdletBinding()]
param(
    [string]$ProjectRoot = 'C:\Users\Admin\Desktop\Scripts Tese 1',
    [string]$Destination = 'C:\Users\Admin\Desktop\Scripts Tese 1\EAAI_REPRODUCIBILITY_PACKAGE_20261015'
)

$ErrorActionPreference = 'Stop'
if (Test-Path -LiteralPath $Destination) {
    throw "Destination already exists: $Destination. Review it before rebuilding; this script never overwrites a package."
}

New-Item -ItemType Directory -Path $Destination | Out-Null
$copied = [System.Collections.Generic.List[object]]::new()

function Copy-PackageFile {
    param([string]$Source, [string]$RelativeDestination)
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) { throw "Required source file is missing: $Source" }
    $target = Join-Path $Destination $RelativeDestination
    $parent = Split-Path -Parent $target
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    Copy-Item -LiteralPath $Source -Destination $target
    $copied.Add([pscustomobject]@{ source = $Source; packaged_path = $RelativeDestination })
}

function Copy-SelectedTree {
    param([string]$SourceDirectory, [string]$RelativeDestination, [string[]]$Extensions)
    if (-not (Test-Path -LiteralPath $SourceDirectory -PathType Container)) { throw "Required source directory is missing: $SourceDirectory" }
    Get-ChildItem -LiteralPath $SourceDirectory -File -Recurse | Where-Object { $Extensions -contains $_.Extension.ToLowerInvariant() } | ForEach-Object {
        $relative = $_.FullName.Substring($SourceDirectory.Length).TrimStart('\')
        Copy-PackageFile -Source $_.FullName -RelativeDestination (Join-Path $RelativeDestination $relative)
    }
}

# Documentation and package instructions.
Copy-PackageFile "$ProjectRoot\Artigo EAAI\reproducibility_README.md" 'README.md'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\reproducibility_REQUIREMENTS.txt" 'requirements.txt'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\reproducibility_RELEASE_CHECKLIST.md" 'PROVENANCE\RELEASE_CHECKLIST.md'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\reproducibility_gitignore" '.gitignore'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\reproducibility_CITATION.cff.template" 'PROVENANCE\CITATION.cff.template'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\EAAI_PRE_SUBMISSION_REVIEW_20261015.md" 'PROVENANCE\EAAI_PRE_SUBMISSION_REVIEW_20261015.md'
Copy-PackageFile "$ProjectRoot\build_eaai_reproducibility_package.ps1" 'PROVENANCE\build_eaai_reproducibility_package.ps1'

# Final data products and provenance.
Copy-SelectedTree "$ProjectRoot\Fase 2\DATASET_CONJUNTO_MCNP63_DEDF_AP_LSDYNA_FINAL_20261002" 'DATA\FINAL_DATASETS' @('.csv', '.json')

# Core analysis source.  Solver executables and original raw output are intentionally excluded.
$codeFiles = @(
    @{ s='Fase 2\integrar_dataset_dedf_ap.py'; d='CODE\PHASE2\integrar_dataset_dedf_ap.py' },
    @{ s='Fase 2\consolidar_dataset_mcnp_lsdyna.py'; d='CODE\PHASE2\consolidar_dataset_mcnp_lsdyna.py' },
    @{ s='Fase 2\postprocess_final_mcnp63_dataset.py'; d='CODE\PHASE2\postprocess_final_mcnp63_dataset.py' },
    @{ s='Fase 2\mcnp63_canonical_validation.py'; d='CODE\PHASE2\mcnp63_canonical_validation.py' },
    @{ s='Fase 3\Fase_3_Metamodelagem_Tese_Final.py'; d='CODE\PHASE3\Fase_3_Metamodelagem_Tese_Final.py' },
    @{ s='Fase 3\feature_engineering_final.py'; d='CODE\PHASE3\feature_engineering_final.py' },
    @{ s='Fase 3\medir_inferencia_candidatos_fase3.py'; d='CODE\PHASE3\medir_inferencia_candidatos_fase3.py' },
    @{ s='Fase 3\tabelas_artigo_fase3.py'; d='CODE\PHASE3\tabelas_artigo_fase3.py' },
    @{ s='Fase 4\Fase_4_Otimizacao_Tese_Final.py'; d='CODE\PHASE4\Fase_4_Otimizacao_Tese_Final.py' },
    @{ s='Fase 4\comparacoes_otimizacao_20261012.py'; d='CODE\PHASE4\comparacoes_otimizacao_20261012.py' },
    @{ s='Fase 4\EXECUTAR_COMPARACOES_20261012.ps1'; d='CODE\PHASE4\EXECUTAR_COMPARACOES_20261012.ps1' },
    @{ s='Fase 5\selecionar_candidatos_validacao_fase5.py'; d='CODE\PHASE5\selecionar_candidatos_validacao_fase5.py' },
    @{ s='Fase 5\preparar_validacao_alta_fidelidade.py'; d='CODE\PHASE5\preparar_validacao_alta_fidelidade.py' },
    @{ s='Fase 5\executar_validacao_mcnp.py'; d='CODE\PHASE5\executar_validacao_mcnp.py' },
    @{ s='Fase 5\consolidar_validacao_alta_fidelidade.py'; d='CODE\PHASE5\consolidar_validacao_alta_fidelidade.py' },
    @{ s='Fase 5\CONVERGENCIA_20261014\preparar_malhas.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\preparar_malhas.py' },
    @{ s='Fase 5\CONVERGENCIA_20261014\executar_malhas.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\executar_malhas.py' },
    @{ s='Fase 5\CONVERGENCIA_20261014\preparar_isolados.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\preparar_isolados.py' },
    @{ s='Fase 5\CONVERGENCIA_20261014\auditar_historicos.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\auditar_historicos.py' },
    @{ s='Fase 5\CONVERGENCIA_20261014\diagnosticar_extremos.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\diagnosticar_extremos.py' },
    @{ s='Fase 5\CONVERGENCIA_9012_ISOLADO\executar_isolados.py'; d='CODE\PHASE5\NUMERICAL_SENSITIVITY\executar_isolados.py' }
)
foreach ($item in $codeFiles) { Copy-PackageFile (Join-Path $ProjectRoot $item.s) $item.d }

# Phase 3 data, frozen splits, predictions, summaries, and operational models.
Copy-SelectedTree "$ProjectRoot\Fase 3\analise_desempenho_modelos (Fase 3) - MCNP63_DEDF_AP" 'DATA\PHASE3' @('.csv', '.json')
Copy-SelectedTree "$ProjectRoot\Fase 3\modelos_treinados (Fase 3) - MCNP63_DEDF_AP" 'MODELS' @('.joblib', '.json')

# Final Phase 4 study and the full equal-budget comparison outputs.  Checkpoints/logs are regenerated and omitted.
Copy-SelectedTree "$ProjectRoot\Fase 4\resultados_otimizacao (Fase 4) - MCNP63_DEDF_AP\fase4_20261003_043316" 'RESULTS\PHASE4\FINAL_OPTIMIZATION' @('.csv', '.json', '.png')
Copy-SelectedTree "$ProjectRoot\Fase 4\comparacoes_otimizacao_20261012\final" 'RESULTS\PHASE4\OPTIMIZATION_COMPARISONS' @('.csv', '.json')

# High-fidelity re-evaluation summaries and numerical-sensitivity diagnostics; never copy raw solver databases.
Copy-SelectedTree "$ProjectRoot\Fase 5\VALIDACAO_ALTA_FIDELIDADE\fase4_20261003_043316_MCNP63_DEDF_AP_30" 'RESULTS\PHASE5\HIGH_FIDELITY_VALIDATION_30' @('.csv', '.json')
Copy-SelectedTree "$ProjectRoot\Fase 5\CONVERGENCIA_20261014" 'RESULTS\PHASE5\NUMERICAL_SENSITIVITY' @('.csv', '.json', '.md')
Copy-SelectedTree "$ProjectRoot\Fase 5\CONVERGENCIA_9012_ISOLADO" 'RESULTS\PHASE5\NUMERICAL_SENSITIVITY_9012' @('.csv', '.json', '.md')

# Submission manuscript and only the figures used by it; build intermediates are excluded.
Copy-PackageFile "$ProjectRoot\Artigo EAAI\manuscript3_EAAI.tex" 'MANUSCRIPT\manuscript3_EAAI.tex'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\manuscript3_EAAI.pdf" 'MANUSCRIPT\manuscript3_EAAI.pdf'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\supplementary_information.tex" 'MANUSCRIPT\supplementary_information.tex'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\supplementary_information.pdf" 'MANUSCRIPT\supplementary_information.pdf'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\references.bib" 'MANUSCRIPT\references.bib'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\elsarticle.cls" 'MANUSCRIPT\elsarticle.cls'
Copy-PackageFile "$ProjectRoot\Artigo EAAI\elsarticle-num.bst" 'MANUSCRIPT\elsarticle-num.bst'
Copy-SelectedTree "$ProjectRoot\Artigo EAAI\figures" 'MANUSCRIPT\figures' @('.pdf', '.png', '.jpg', '.jpeg', '.eps')
Copy-SelectedTree "$ProjectRoot\Artigo EAAI\revision_20261012" 'MANUSCRIPT\revision_20261012' @('.pdf', '.png', '.jpg', '.jpeg', '.eps')
Copy-SelectedTree "$ProjectRoot\Artigo EAAI\revision_20261014" 'MANUSCRIPT\revision_20261014' @('.pdf', '.png', '.jpg', '.jpeg', '.eps')

# Final immutable inventory and hashes.
$records = Get-ChildItem -LiteralPath $Destination -File -Recurse | Sort-Object FullName | ForEach-Object {
    $relative = $_.FullName.Substring($Destination.Length).TrimStart('\')
    $hash = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash
    [pscustomobject]@{ packaged_path = $relative; bytes = $_.Length; sha256 = $hash }
}
$records | Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $Destination 'PROVENANCE\sha256_manifest.csv')
$copied | Export-Csv -NoTypeInformation -Encoding utf8 (Join-Path $Destination 'PROVENANCE\source_to_package_map.csv')

[pscustomobject]@{
    package_root = $Destination
    files_hashed = $records.Count
    bytes = ($records | Measure-Object -Property bytes -Sum).Sum
    gigabytes = [math]::Round((($records | Measure-Object -Property bytes -Sum).Sum / 1GB), 3)
} | ConvertTo-Json | Set-Content -Encoding utf8 (Join-Path $Destination 'PROVENANCE\package_build_summary.json')

Write-Host "Package created: $Destination"
Write-Host "Files hashed: $($records.Count)"
