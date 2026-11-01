# Publish the Dual-Lock repository to GitHub (run in PowerShell; requires git and the GitHub CLI `gh`, already logged in as sandlerleon)
$repoDir = "C:\Users\Leon\Downloads\DualLock_work\repo"
$slug    = "sandlerleon/dual-lock-oral-peptide-mc"
Set-Location $repoDir

function Check($what) { if ($LASTEXITCODE -ne 0) { Write-Host "FAILED: $what" -ForegroundColor Red; exit 1 } }

# 1. create the public repository if it does not exist yet (a 'not found' message from gh is expected the first time)
$exists = $false
$prev = $ErrorActionPreference; $ErrorActionPreference = "Continue"
$null = & gh repo view $slug 2>&1
if ($LASTEXITCODE -eq 0) { $exists = $true }
$ErrorActionPreference = $prev
if (-not $exists) {
    gh repo create $slug --public --description "Dual-Lock oral-peptide Monte Carlo: mechanistic design-space model, tests and results (Int. J. Pharm. manuscript)"
    Check "gh repo create"
}

# 2. remote, push branch and the v1.0.0 tag
$ErrorActionPreference = "Continue"
git remote remove origin 2>$null
$ErrorActionPreference = $prev
git remote add origin "https://github.com/$slug.git"; Check "git remote add"
git branch -M main
git push -u origin main; Check "git push main"
git push origin v1.0.0 --force; Check "git push tag"

# 3. release with the manuscript PDFs
gh release create v1.0.0 `
    "manuscript\DualLock_IJP_Manuscript.pdf" `
    "manuscript\DualLock_IJP_Supplementary_Information.pdf" `
    "manuscript\Cover_Letter_IJP.pdf" `
    --repo $slug --title "v1.0.0 - Dual-Lock oral-peptide Monte Carlo" `
    --notes "Code, tests, results and manuscript for the International Journal of Pharmaceutics submission. Zenodo: software 10.5281/zenodo.23864330, preprint 10.5281/zenodo.23864327."
Check "gh release create"

Write-Host "Done: https://github.com/$slug" -ForegroundColor Green
