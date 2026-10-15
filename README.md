# Sensitivity analysis for a meta-analysis of APIM paths

How large do pooled Actor–Partner Interdependence Mediation Model (APIMeM)
paths have to be before they become statistically significant, given the
studies actually available? This repository holds the simulation that answers
that question, for both an uncorrected meta-analysis and one corrected for
the unreliability of the questionnaire measures.

Nothing is simulated at the dyad level. Each study enters as its 12
standardized path estimates and their full 12 × 12 sampling covariance
matrix, so the simulation reproduces the precision of the real design without
reusing the raw dyadic data.

## What the script does

`apim_meta_sensitivity.R` runs the whole analysis end to end:

1. Reads study-level estimates and sampling covariances from
   `APIM_meta_analysis_results.rds` (or from a table typed into the script —
   see `CFG$input_mode`).
2. Simulates pooled estimates under four levels of between-study
   heterogeneity (I² = 0, 25, 50, 75%).
3. Fits six models to every simulated dataset:
   - `metafor::rma()` — common effect (z), random effects (z), Knapp–Hartung,
     Knapp–Hartung ad hoc
   - `metafor::rma.mv()` — multilevel (partner-specific estimates nested in
     studies), with model-based z and cluster-robust CR2 intervals from
     `clubSandwich`
4. Derives power curves, minimum detectable effects, critical values, Type I
   error rates, and joint-significance power surfaces for the indirect
   effects.
5. Writes the figures, summary tables and an APA-formatted Word report.

Because every model is an intercept model with a fixed sampling covariance
matrix, the standard errors and critical values do not depend on the true
effect. Each replication is therefore simulated once under a null effect and
reused across the whole grid, which is what makes complete power surfaces
affordable. Section 6 of the script verifies that invariance numerically.

## Running it

```r
source("apim_meta_sensitivity.R")
```

R ≥ 4.2 with `metafor`, `clubSandwich`, `ggplot2`, `flextable`, `officer`,
`MASS` and `parallel`. At the default 1,000 replications it takes roughly ten
minutes on one core.

Settings live in the `CFG` list at the top of the script. To change them
without editing the file, define `CFG_OVERRIDE` first:

```r
CFG_OVERRIDE <- list(R = 300, het_values = c(0, .5), out_prefix = "quick")
source("apim_meta_sensitivity.R")
```

`CFG$k_replicates = c(1, 2, 5)` adds hypothetical studies with the same
sampling covariances as the observed ones and produces an extra figure
showing how the minimum detectable effect falls as studies accumulate.

## Files

| File | What it is |
|---|---|
| `apim_meta_sensitivity.R` | the analysis script |
| `APIM_meta_analysis_results.rds` | study-level input: estimates and sampling covariances |
| `APIM_sensitivity_report.docx` | report: four tables and seven figures |
| `APIM_sensitivity_fig1–2_*.png` | power curves for the direct paths |
| `APIM_sensitivity_fig3_mde_direct.png` | minimum detectable effect by heterogeneity and method |
| `APIM_sensitivity_fig4_critical_estimate.png` | smallest pooled estimate that would reach significance |
| `APIM_sensitivity_fig5_type1_error.png` | empirical Type I error rates |
| `APIM_sensitivity_fig6–7_indirect_*.png` | 80% power contours over the (a, b) grid |
| `APIM_sensitivity_direct_summary.csv` | minimum detectable effects, critical values, Type I error |
| `APIM_sensitivity_indirect_summary.csv` | smallest detectable indirect effects |
| `APIM_sensitivity_observed_critical.csv` | observed pooled estimates and their critical values |
| `APIM_sensitivity_mc_vs_joint.csv` | joint significance vs Monte Carlo CI check |

Two generated files are not tracked, because the script recreates them and
they are large: the full power-curve grid
(`APIM_sensitivity_direct_power_curves.csv`, about 10 MB) and the raw results
object (`APIM_sensitivity_results.rds`).

The input `.rds` is produced by the companion script `apim_meta_analysis.R`,
which fits the APIMeM to each study and pools the paths. It is not part of
this repository.

## Reading the results

With only two studies the methods differ sharply, and the differences are not
cosmetic:

- The z-based tests (common effect, random effects, multilevel) detect pooled
  effects of roughly .11 and above when heterogeneity is absent, rising to
  about .22 at I² = 75%; the reliability-corrected paths need about .15 to .34
  on their own scale.
- Knapp–Hartung and CR2 test on about one degree of freedom. Their critical
  values swing with how closely the two studies happen to agree, and the
  ad hoc Knapp–Hartung variant never rejects at all in this design. Treat the
  robust results as descriptive until more studies are available.
- Under heterogeneity the common-effect test rejects well above the nominal
  rate, because it is testing the weighted average of these particular
  studies rather than a population mean.

## License

MIT (see `LICENSE`). The copyright line reads "The Authors" so the repository
can be shared during masked review; replace it with your name afterwards.
