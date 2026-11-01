# Publish the Dual-Lock repository to GitHub (run in PowerShell; requires git and the GitHub CLI `gh`, already logged in as sandlerleon)
$ErrorActionPreference = "Stop"
$repoDir = "C:\Users\Leon\Downloads\DualLock_work\repo"
$slug    = "sandlerleon/dual-lock-oral-peptide-mc"
Set-Location $repoDir

# 1. create the public repository (skipped if it already exists)
gh repo view $slug 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    gh repo create $slug --public --description "Dual-Lock oral-peptide Monte Carlo: mechanistic design-space model, tests and results (Int. J. Pharm. manuscript)"
}

# 2. remote, push branch and the v1.0.0 tag
git remote remove origin 2>$null
git remote add origin "https://github.com/$slug.git"
git branch -M main
git push -u origin main
git push origin v1.0.0 --force

# 3. release with the manuscript PDFs
gh release create v1.0.0 `
    "manuscript\DualLock_IJP_Manuscript.pdf" `
    "manuscript\DualLock_IJP_Supplementary_Information.pdf" `
    "manuscript\Cover_Letter_IJP.pdf" `
    --repo $slug --title "v1.0.0 - Dual-Lock oral-peptide Monte Carlo" `
    --notes "Code, tests, results and manuscript for the International Journal of Pharmaceutics submission. Zenodo: software 10.5281/zenodo.23864330, preprint 10.5281/zenodo.23864327."

Write-Host "Done: https://github.com/$slug"
