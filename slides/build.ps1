$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    $build = Join-Path $env:TEMP ('saw_beamer_' + [guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $build | Out-Null
    1..2 | ForEach-Object {
        pdflatex -interaction=nonstopmode -halt-on-error "-output-directory=$build" tesi_saw.tex
        if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed. Log: $build" }
    }
    Copy-Item -LiteralPath (Join-Path $build 'tesi_saw.pdf') -Destination 'tesi_saw.pdf' -Force
    Write-Output "PDF: $PSScriptRoot\tesi_saw.pdf"
    Write-Output "Build log: $build"
} finally { Pop-Location }
