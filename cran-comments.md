## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission.

Changes in this version include:

* Replaced `method.tau` with separate arguments for HK, RE and heterogeneity estimation in `confMeta()`.
* Improved handling of `sm` in the internal `meta::metagen()` call.
* Added `ref_labels` to `autoplot.confMeta()`.
* Improved forest plot y-axis spacing.