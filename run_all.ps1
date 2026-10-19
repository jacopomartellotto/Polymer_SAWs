param(
    [int]$Jobs = 8,
    [switch]$Overwrite,
    [switch]$DryRun,
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$arguments = @("-B", "scripts/run_pipeline.py", "--jobs", $Jobs)
if ($Overwrite) { $arguments += "--overwrite" }
if ($DryRun) { $arguments += "--dry-run" }

& $Python @arguments
if ($LASTEXITCODE -ne 0) {
    throw "Pipeline failed with exit code $LASTEXITCODE"
}
