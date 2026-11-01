# Dual-Lock oral-peptide Monte Carlo

Leon Sandler, Independent Researcher (Northbrook, Illinois, USA) — sandler.leon@gmail.com — ORCID 0009-0007-4584-808X

Code, tests, results and manuscript for *"Mechanistic Monte Carlo assessment of co-localized permeation enhancement and efflux modulation for oral peptide delivery,"*
prepared for the **International Journal of Pharmaceutics** (Elsevier).

Archived on Zenodo: software doi:[10.5281/zenodo.23864330](https://doi.org/10.5281/zenodo.23864330) (all versions; version 1.0.0 is 10.5281/zenodo.23864331), manuscript preprint doi:[10.5281/zenodo.23864327](https://doi.org/10.5281/zenodo.23864327) (all versions; v1 is 10.5281/zenodo.23864328).

**This is an in-silico design-space study.** It contains no experimental data. Every probability is conditional on author-assigned, evidence-graded priors and on the model structure;
none is an empirical probability of success. The reference arm is calibrated to a published baseline, so it is not a validation.

## What it does

A lumped mechanistic model of the micro-niche between a tablet and the mucosa (peptide, permeation enhancer and efflux modulator pools; Hill enhancement of permeability; efflux ratio reduced by the
modulator; proteolysis, wash-out and enhancer clearance; damage-hours as an exposure metric) is propagated through evidence-graded priors. Beyond the first version of the model (`oral_peptide_mc.py`, unchanged apart from a tie-correct rank correlation) the extension adds:

* common-random-number (paired) arms, ablations and negative controls; Gaussian-copula correlated priors;
* convergence analysis (running estimates, batch means, 20 replicate runs, time-step and horizon studies);
* global variance-based sensitivity (Saltelli / Jansen estimators, single inputs and groups) and uncertainty decomposition;
* benefit-risk (Pareto) dose sweeps and falsification-region maps;
* structural uncertainty: a stagnant-layer variant and a spatial diffusion-reaction model (implicit finite volumes), each recalibrated to the same baseline;
* observation-conditioned analysis (rejection ABC) and comparison with published human bioavailability data;
* revision analyses: uptake-constrained model (stagnant layer + peptide-solubility limit), decomposition of the matrix-only gain, biologically constrained recalibration, reversibility-endpoint scenarios (C_tox/EC50, damage-hours and recovery limits) and out-of-sample evaluation of an optimized design (train/test split).

## Layout

| Path | Contents |
|---|---|
| `code/oral_peptide_mc.py` | the original lumped model, priors, arms and 9 self-tests (`--test`, `--run`) |
| `code/dualock_lib.py` | uniform-variate sampler, copula, stagnant-layer and spatial models, Sobol estimators |
| `code/run_extended.py` | the analysis stages (`equiv paired converge corr sobol pareto fals struct abc valid constrained revers optim decomp`) and 19 verification tests (`--test`) |
| `code/make_figures.py`, `code/make_tables.py` | figures and Supplementary CSV tables from the result files |
| `tests/test_all.py` | pytest wrapper of both test suites |
| `results/` | `results_v1.json` (original run), `ext_*.json` (all numerical results), `tables/` (CSV), `run_ext.log` |
| `figures/` | the figures of the manuscript and supplement |
| `manuscript/` | manuscript, cover letter, highlights and Supplementary Information (docx, pdf) with their build scripts in `manuscript/build/` |

## Reproduce

```bash
pip install -r requirements.txt
cd code
python oral_peptide_mc.py --test          # 9 self-tests of the original model (about 25 s)
python run_extended.py --test             # 19 verification tests of the extension (about 1 min)
python oral_peptide_mc.py --run           # original analysis -> results/results.json (about 2 min); headline hash 4df279944ef2355e
python run_extended.py                    # all extended stages (about 30 min on one core); or name stages: python run_extended.py sobol fals
python make_figures.py results_ext        # figures
python make_tables.py results_ext         # Supplementary CSV tables
```

Seeds are fixed (original 20261031; extension 20261101 + offsets), so a re-run regenerates every number. Reported with Python 3.11.9, NumPy 2.4.6, SciPy 1.17.1 and Matplotlib 3.11.1.

## Reproduction of the first version

The model of the first version was re-run unchanged in an empty directory: the headline hash (4df279944ef2355e) and all arm results were identical; 18 of 8,683 numbers differed, all Spearman
coefficients of the sensitivity block (by at most 0.03), because tied ranks (a binary success indicator; the efflux ratio of non-substrates) were broken arbitrarily. The tie-correct version is used here.

## Citation

See `CITATION.cff`. MIT licence.
