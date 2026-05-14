# confMeta 0.1.1

* Replaced `method.tau` with three separate `method.tau` arguments to independently control tau estimation for HK, RE and heterogeneity estimation in `confMeta()`.
* In the internal `meta::metagen()` call used by `confMeta()`, replaced `sm = "MD"` with `sm = ""` to improve clarity of the code.
* Added `ref_labels` to `autoplot.confMeta()` to allow custom display labels.
* Improved forest plot by increasing the y-axis upper limit by half a polygon height.

# confMeta 0.1.0

* Initial CRAN submission.