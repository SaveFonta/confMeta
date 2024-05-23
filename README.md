# confMeta

The `confMeta` package implements methods related to meta-analysis. The
functions that `confMeta` provides, can be categorized roughly into three
categories:

* An S3 class `confMeta` and a plot method for it
* *p*-value functions
* Heterogeneity adjustments for standard errors

## Installation

Currently, the package is developed on
[Github](https://github.com/felix-hof/confMeta). Thus, the easiest way to
install it is via the `remotes` package. If the `remotes` package is not
installed, you can run the following line of code to do so.

```r
install.packages("remotes")
```

Once `remotes` has been installed, you can install `confMeta` by running

```r
remotes::install_github("felix-hof/confMeta")
```

Once installed, the package can be used by loading it into the `R` session. This
can be achieved by running

```r
library(confMeta)
```

## Usage examples


### Simulating individual studies
Consider the hypothetical scenario where we have *n = 3* individual studies
that should be combined into a single confidence interval using a confidence
level of *1 - alpha = 0.95*. We can simulate these by running the following code

```r
n <- 3
conf_level <- 0.95
estimates <- rnorm(n)
SEs <- rgamma(n, 5, 5)
```

Here, the object `estimates` contains the individual study estimates, whereas
the object `SEs` contains the corresponding standard errors.

### Creating the confMeta object

With these individual studies, a `confMeta` object can be created. However, this
requires the specification of a *p*-value function, i.e. a method, that takes
the individual studies (argument `estimates`), their standard errors
(argument `SEs`), and the mean under the null-hypothesis (argument `mu`)
as input and returns the corresponding *p*-value at the specified mean value.
The `confMeta` package provides implementations for the following *p*-value
functions

* Harmonic mean (Function: `p_hmean`)
* Wilkinson (Function: `p_wilkinson`)
* Pearson (Function: `p_pearson`)
* Edgington (Function: `p_edgington`)
* Fisher (Function: `p_fisher`)
* Tippett (Function: `p_tippett`)

In this example, we choose Edgington's method. Thus, we can create the 
`confMeta` object as follows

```r
cm <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    fun = p_edgington,
    fun_name = "Edgington"
)
```

As the variable `cm` now contains the `confMeta` object, we can inspect it by
running the following code 

```r
# See what elements it has
names(cm)

# Check out the combined confidence interval(s)
cm$joint_cis
```

### Visualizations

The package also contains an `autoplot` method that can be used to visualize the
*p*-value function. The documentation for this function can be inspected by
running

```r
?autoplot.confMeta
```

The method provides essentially two plots, one showing the *p*-value
function and one constructing a forest plot. Which one is returned can be
specified using the `type` argument.

```r
# show the p-value function
autoplot(cm, type = "p")

# show the forest plot
autoplot(cm, type = "forest")

# show both
autoplot(cm, type = c("p", "forest"))
```

You can also compare different *p*-value functions with each other. In order to
illustrate how this works, we need a second `confMeta` object, that uses a
different *p*-value function.

```r
cm2 <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    fun = p_fisher,
    fun_name = "Fisher"
)
```

Now, we can compare the two *p*-value functions to each other with in the
following way

```r
# show the p-value function
autoplot(cm, cm2, type = "p")

# show the forest plot
autoplot(cm, cm2, type = "forest")

# show both
autoplot(cm, cm2, type = c("p", "forest"))
```
