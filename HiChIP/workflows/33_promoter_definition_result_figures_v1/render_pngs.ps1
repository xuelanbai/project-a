param(
    [Parameter(Mandatory=$true)][string]$OutputDirectory,
    [int]$Resolution = 220
)

$ErrorActionPreference = 'Stop'
$pdftoppm = 'C:\Users\Luca\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe'
if (-not (Test-Path -LiteralPath $pdftoppm)) {
    throw "Bundled pdftoppm was not found: $pdftoppm"
}
Get-ChildItem -LiteralPath $OutputDirectory -Filter '*.pdf' | ForEach-Object {
    $prefix = Join-Path $OutputDirectory $_.BaseName
    & $pdftoppm -png -singlefile -r $Resolution $_.FullName $prefix
    if ($LASTEXITCODE -ne 0) {
        throw "PDF rendering failed: $($_.FullName)"
    }
}
