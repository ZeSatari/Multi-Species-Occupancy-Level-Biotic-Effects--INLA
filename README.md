# Spatio-temporal occupancy models with observation-level biotic information

Analysis code for *Observation-Level Biotic Information in Single- and
Multi-Species Spatio-Temporal Occupancy Models with INLA-SPDE*.

Every numerical claim in the manuscript is checked against the analysis output
by `R/90_verify_claims.R`. Run it after any rerun; it reports which statement
in the text is affected if a number changes.

## Environment

| | |
|---|---|
| R | 4.5.3 |
| R-INLA | 25.10.19 (from `https://inla.r-inla-download.org/R/stable`) |
| Other packages | `fmesher`, `sf`, `terra`, `dplyr`, `MASS`, `future.apply`, `sn`, `scoringRules`, `spOccupancy`, `ggplot2`, `patchwork` |
| Platform | Windows 11, 4 workers |

`sn` is required by `inla.posterior.sample()` and `inla.hyperpar.sample()`;
without it the fits run but no posterior draws can be taken. `scoringRules` is
used only to cross-check the CRPS implementation.

**Fix one INLA version and reproduce every number under it.** The `occupancy`
likelihood family is recent enough that behaviour may differ between releases,
and results from two versions must not be mixed in one table.
`record_session()` writes the version alongside each set of results.

## Run order

```
R/00_setup.R            packages, priors, options, shared helpers
R/10_metrics.R          bias, RMSE, coverage, CRPS
R/20_simulate_single.R  single-species generator, four regimes
R/30_fit_occupancy.R    the single fitting entry point

R/70_strategy_check.R     approximation adequacy      -> results/strategy/
R/40_run_single_sim.R     single-species simulation   -> results/single/
R/50_run_msom_sim_v2.R    multi-species simulation    -> results/msom/
R/61_stability_check.R    empirical fit reproducibility
R/60_real_data.R          Hubbard Brook application   -> results/real/
R/80_figure1.R            Figure 1                    -> figures/
R/90_verify_claims.R      checks every numerical claim in the paper
```

Diagnostics and checks, kept because each is cited in the paper:

```
R/62_missing_sensitivity.R  treatment of absent co-occurring records
R/63_avg_occupancy.R        average occupancy vs observed frequency
R/64_mesh_adequacy.R        three mesh resolutions compared
R/65_site_cv.R              site-held-out cross-validation
R/66_prior_sensitivity.R    species-effect prior, four specifications
R/67_sigma_gamma_check.R    reproducibility of sigma_gamma
R/68_bias_source.R          slope bias, likelihood in isolation
R/69_bias_field.R           slope bias, field and covariate form
R/99_diag_group.R           the `A.local` grouping convention
R/99_diag_confounding.R     covariate/field separability
```

`R/68_bias_source.R` also contains the only comparison in this work against a
method with no approximation: at perfect detection the occupancy model reduces
to logistic regression, and R-INLA agrees with `glm()` to within 0.002 on
identical data. That covers the mixture likelihood, not the SPDE representation.

`R/50_run_msom_sim.R` and `R/40_run_single_sim-old.R` are superseded and
should be deleted.
`R/50_run_msom_sim_v3.R` is the homogeneous-slope diagnostic reported in
Appendix S1. `R/50_run_msom_sim_v4.R` was written but not used; see
"Not addressed" below.

Set `RUN_FULL <- FALSE` before sourcing `40_` or `50_v2` to load the
configuration and replicate function without launching the full study — useful
for timing one replicate first.

## Runtimes

Measured on the platform above.

| Script | Fits | Time |
|---|---|---|
| `70_strategy_check` | 3 | ~45 min |
| `40_run_single_sim` | 1600 | ~10 h (4 workers) |
| `50_run_msom_sim_v2` | 1300 | ~11 h (4 workers) |
| `60_real_data` | 16 | ~1 h |

`50_v2` writes a checkpoint every ten replicates to
`results/msom/checkpoints/`. If interrupted, rerun the same command and it
resumes from the last completed chunk; seeds are `SEED_BASE + r`, so a
replicate run later is identical to the one it replaces. Delete the checkpoint
directory to force a clean run.

Each replicate of `50_v2` holds ~2.6 GB. Do not use `detectCores()` on a
hybrid P-core/E-core CPU: workers on efficiency cores run several times slower
and `future_lapply` waits for the slowest. Set `N_WORKERS` by hand.

The runs reported in the paper used four workers on an Intel Core Ultra 5 225U
with 16 GB. `detectCores()` reports 12 physical and 14 logical cores there, but
only two of those carry hyperthreads; the rest are efficiency cores. Four was
also close to the memory ceiling at 2.6 GB per replicate. Both constraints
pointed to the same number, which is why it is hard-coded rather than
detected.

## Findings that shape the code

**Mesh resolution must be chosen relative to the information the data can
support, not the geometry of the study area.** An initial mesh of 1724 nodes
gave a latent field of 15,516 elements against 3,357 observation rows.
Repeated fits of the identical model then returned intercepts spanning 0.183
on the logit scale, run times varying between 1,900 and 13,100 seconds, and
convergence failures in up to 80% of attempts. Coarsening to 283 nodes — fewer
field elements than observations — reduced the spread to 0.001 with uniform
convergence. `60_real_data.R` warns if the ratio exceeds 2. Because the
occupancy log-likelihood is not log-concave, verify that repeated fits
reproduce rather than assuming they do; `61_stability_check.R` does this.

**The approximation was checked, not assumed.** Refitting under Gaussian,
simplified Laplace and full Laplace strategies gave occupancy posterior means
agreeing to three decimal places, with a largest discrepancy of 0.042 posterior
standard deviations. This is an internal check, not an external benchmark: no
comparison against MCMC or `spOccupancy` was made.

**`int.strategy = "eb"` is not used.** Empirical Bayes fixes the
hyperparameters at their posterior mode, so credible intervals understate
uncertainty and coverage is not interpretable. All reported results use
`"ccd"`.

**A quadratic covariate attenuates its own coefficient.** With `x ~ N(0,1)`,
`x^2` is chi-squared with one degree of freedom; its right tail drives psi
toward one, where the data carry little information about the slope.
Diagnostics with no spatial structure at n = 5000 recovered 0.697 for a linear
covariate against a truth of 0.700, but 0.456 for raw `x^2` and 0.614 after
standardisation. The multi-species simulation therefore uses a linear
covariate so that residual bias is attributable to the hierarchical structure.
The empirical application retains a quadratic elevation term.

**The "unused groups" warning is a false positive with `A.local`.** INLA counts
groups from the index vector, which is all `NA` by design. Verified against a
simulated AR(1) field: the correct specification recovers rho = 0.794 for a
generating value of 0.80, while the specification that suppresses the warning
returns 0.000 and a substantially worse fit (`R/99_diag_group.R`). Do not
"fix" it.

**The detection model is capped at ten coefficients**, including the
intercept, because the occupancy family declares a fixed number of
hyperparameters (`theta1` ... `theta10`). The cap belongs to the current
implementation, not to the INLA-SPDE approach, and could be raised. Raising it
would not remove the underlying constraint: detection coefficients enter as
hyperparameters rather than as latent field elements, and numerical
integration over them becomes computationally infeasible when their dimension
is large (Belmont et al. 2024), so the approximation deteriorates well before
the cap is reached. `build_det_X()` errors at the cap and warns within two of
it.

## Defects corrected from the earlier pipeline

| Old behaviour | Consequence |
|---|---|
| `mvrnorm(mu = rep(-0.5, 12))` plus `mu_beta0 = -0.5` | true community intercept was −1.0, reported as −0.5 |
| Generated corr +0.5, fitted `Cmatrix = I + J` (corr −1/12) | generation and fitting misspecified in opposite directions |
| `scale_x_s^2` in a formula | `^` is the crossing operator; the quadratic term was never fitted |
| `as.vector(x_covariate)` assigned by position | lengths matched, ordering did not |
| `inla.group.cv(result = model_xs, ...)` | undefined object; a published ULOOCV value came from nowhere |
| CAWA block never subset the species | one table row came from a model that could not run |
| No `constr` on `f(species_id, "iid")` | beta0 and gamma_i not separately identified; CrI 5.3x too wide |
| WAIC/n compared across different response sets | the column tracked prevalence, not model quality |
| "MSE" = squared posterior SD of one fit | measures posterior spread, not distance from truth |
| R = 1 | bias, coverage and RMSE are undefined |

Two further defects were introduced during the rewrite and caught by testing:
`species_intercepts()` indexed posterior samples by position rather than name,
attaching species effects to the wrong species; and a `constr = TRUE` random
slope without a companion fixed effect forced the community slope to zero.
Both produced plausible output and neither raised an error.

## Verification

`R/90_verify_claims.R` restates each numerical claim in the manuscript as a
testable expression and evaluates it against the saved results. It prints a
PASS/FAIL table and stops on any failure.

Adding a number to the manuscript means adding a `claim()` line. Two errors
were caught this way, both in text drafted from a single replicate rather than
from the full run.

```r
source("R/90_verify_claims.R")
# 56 of 56 claims verified.
```

## Outputs

```
results/strategy/   approximation comparison, timing, marginal overlays
results/single/     per-replicate parameters and fit statistics, Tables S1-S3
results/msom/       community and species-level recovery, Tables S4-S5
results/real/       detection coefficients, species table, confounding
                    diagnostic, stability check, fitted models
figures/            figure1.pdf, figure1.png
results/claim_verification.csv
```

Each results directory carries `sessionInfo.txt` and
`priors_and_options.rds`.

## Data

The empirical data are the `hbefTrends` object in the `spOccupancy` package,
from the Hubbard Brook Ecosystem Study. No data files are redistributed here.

## Not addressed

- **No MCMC benchmark.** The occupancy log-likelihood is not log-concave, so
  the accuracy of the approximation for this model class is established only
  by the internal strategy comparison. Validation against NIMBLE or Stan
  remains outstanding.
- **No comparison with `spOccupancy`**, from which the empirical data are
  drawn. No computational advantage is claimed in its absence.
- **Prior sensitivity.** `PRIORS_SENS` is defined in `00_setup.R` but no
  driver uses it.
- **Species-specific loadings** on the shared field are discussed in the paper
  as an extension but not implemented.
- **Source of the slope bias, not found.** `R/68_bias_source.R` and
  `R/69_bias_field.R` were written to locate the negative bias in the
  environmental slope (-0.024 single-species, -0.050 multi-species). Both ran
  to completion; neither reproduced it. Stripped to one species with a linear
  covariate, no field and no random effects, the estimator is unbiased across
  sample size (250-5000), detection probability (0.15-1.0), visits (2-5) and
  all three approximation strategies, with coverage near nominal; at p = 1,
  R-INLA and `glm()` agree to within 0.002. A quadratic covariate produces a
  bias of -0.21, an order of magnitude too large. Fitting the latent field
  reduces the bias rather than causing it. What remains is unexplained and the
  manuscript says so.

  One incidental finding: at K = 10 visits every fit failed (`n_ok = 0` in
  `results/bias/bias_source.csv`). The detection design has K columns per
  covariate and something in that path does not scale; not investigated.

- **Heterogeneous-prevalence simulation, attempted and abandoned.**
  `R/50_run_msom_sim_v4.R` is a complete driver for the configuration the
  empirical data exhibit: species effects drawn to reproduce the observed
  spread, with occupancy and detectability correlated at 0.9 as they are at
  Hubbard Brook (the observed correlation on the logit scale is 0.997). It was
  never run past a single test replicate.

  On that replicate the hierarchical estimates for the sparsest species did not
  shrink toward the community mean but overshot the generating values in the
  opposite direction: a species generated at -3.67 was estimated at -5.65,
  while the independent fit gave -0.39. The two estimators erred in opposite
  directions and neither was close. Widening the prior on the species-effect
  standard deviation, which the generating spread (sd 1.38) placed far in the
  tail of `PRIORS$sigma_gamma` (P(sigma > 1) = 0.01), changed nothing; the
  widened version is retained in the script as `PRIORS_V4_GAMMA`.

  The full run was not attempted because comparing error metrics between two
  estimators, one behaving in a way we could not explain, would not have been
  informative. Whether the behaviour is specific to this generator or says
  something about the model in that regime is unresolved, and settling it is a
  prerequisite for the run rather than a detail to sort out afterwards. The
  manuscript reports the attempt in its Limitations.
