# GeneticScreen

Pre- and post-processing for the analysis of high-throughput genetic screens using matrix linear models. See the associated paper, ["Matrix linear models for high-throughput chemical genetic screens"](http://dx.doi.org/10.1534/genetics.119.302299), for more details. S-scores implemented based on Collins et al. (2006) <sup>[1](#myfootnote1)</sup>. 

`GeneticScreen` is an extension of the [`matrixLM`](https://github.com/janewliang/matrixLM.jl) package, which provides core functions for closed-form least squares estimates for matrix linear models. 

## Installation

The `GeneticScreen` package can be installed by running: 

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/janewliang/GeneticScreen.jl", rev="master"))
```

`GeneticScreen` was developed in [Julia v1.1](https://julialang.org/downloads/). 

## Usage

```
using GeneticScreen
```

The `GeneticScreen` package extends [`matrixLM`](https://github.com/janewliang/matrixLM.jl), so all functionality of the `matrixLM` functions is preserved. The README for `matrixLM` provides examples of usage for the core least squares estimation functions that may be of interest to the user. 

In this illustrative example, consider simulated data from a tiny genetic screening experiment run on a single 4x6 plate. There were 5 experimental conditions (A-E), each with 4 monotonically ordered dosage concentration levels (1-4) replicated 2 times. There were 8 mutants, each replicated 3 times per plate. The simulated data is stored in the "data/" directory as CSV files. 

```
using DataFrames
using CSV
using Random

# Generate 5 conditions (A-E), each with 4 ordered dosage concentration 
# levels (1-4)
Xdos_df = repeat(DataFrame(cond=repeat(["A", "B", "C", "D", "E"], inner=4), 
                           conc=repeat(1:4, 5)), 2)
# Write Xdos_df to CSV
CSV.write("./data/Xdos_df.csv", Xdos_df)

# Create another DataFrame that concatenates the conditions and concentrations
X_df = repeat(DataFrame(cond_conc=[string(Xdos_df[:cond][i], Xdos_df[:conc][i]) 
                        for i in 1:20]), 2)
# Write X_df to CSV
CSV.write("./data/X_df.csv", X_df)

# Generate 8 mutants, each replicated 3 times
Z_df = DataFrame(mut=repeat(1:8, 3))
# Write Z_df to CSV
CSV.write("./data/Z_df.csv", Z_df)

# Generate the same 8 mutants with 3 replicates, and also include the spatial 
# row and column positions on the 4x6 plate
Zspat_df = DataFrame(mut=repeat(1:8, 3), 
                     row=repeat(1:4, inner=6), col=repeat(1:6, 4))
# Write Zspat_df to CSV
CSV.write("./data/Zspat_df.csv", Zspat_df)

# Use `contr` to get sum contrasts for the concatenated 
# condition-concentrations in X_df and the mutants in Z_df, and convert them 
# into 2d arrays
X = convert(Array{Float64,2}, contr(X_df, [:cond_conc], ["sum"]))
Z = convert(Array{Float64,2}, contr(Z_df, [:mut], ["sum"]))

# Total number of condition replicates (rows of Y)
n = size(X)[1]
# Total number of mutant replicates per plate (columns of Y)
m = size(Z)[1]
# Number of conditions
p = size(X)[2]
# Number of mutants
q = size(Z)[2]

# Randomly generate response variable Y
Random.seed!(1)
B = rand(-5:5,p,q)
E = randn(n,m)
Y = X*B*transpose(Z)+E
# Write Y to CSV
CSV.write("./data/Y.csv", DataFrame(Y))
```

The `read_plate` function is designed to simplify the construction of a `RawData` object for genetic screening data. (The `RawData` is needed to obtain least squares estimates for matrix linear models.) Users can specify the paths to flat files where their X (experimental conditions and row covariates), Y (colony size/response), and Z (mutants and column covariates) matrices are stored. `read_plate` will read in the data using `CSV.read`; additional keyword arguments to be passed into `CSV.read` for all three files can be specified. By default, `CSV.read` assumes that the first row is a header and that the delimiter is ','. 

If one of the columns in the X matrix is a categorical variable for the experimental conditions, `read_plate` can convert it to sum or treatment contrasts using `contr` (see the `matrixLM` package for more information on `contr`). The condition variable name should be specified using the `XCVar` argument and the type of contrast should be specified using `XCType`. If one of the columns in the Z matrix is a categorical variable for the mutants, the analogous argument names are `ZCVar` and `ZCType`. Additional variables in X and Z are retained in the resulting `RawData` object. In the example below, we use an X matrix that concatenates the conditions with their concentrations, thereby creating 24 distinct "condition-concentrations" (i.e. condition A at concentration 1 is treated as a different "condition" from condition A at concentration 2). By default, both the condition-concentrations in X and the mutants in Z are coded as treatment contrasts, but they are coded as sum contrasts in this example. 

```
dat1 = read_plate("./data/X_df.csv", "./data/Y.csv", "./data/Z_df.csv", 
                  XCVar=:cond_conc, XCType="sum", 
                  ZCVar=:mut, ZCType="sum")
```

If monotonicity is a reasonable assumption, it may be of interest to consider the concentrations of the different conditions as dosage slopes. If `isXDos=true`, `read_plate` will code the concentrations in `XConcentrationVar` as dosage slopes for each condition in `XConditionVar`. This is done through a call to the `get_dose_slopes` function, which assumes that concentrations are ordered within dosage slopes. 

```
dat2 = read_plate("./data/Xdos_df.csv", "./data/Y.csv", "./data/Z_df.csv", 
                  isXDos=true, XConditionVar=:cond, XConcentrationVar=:conc, 
                  ZCVar=:mut, ZCType="sum")
```

Colonies grown on the edge of the plate tend to be larger because they do not compete for as many resources. One way to try to account for these spatial differences is to include power and cross-product terms for the (centered) row and column indices. In this example, the rows are indexed from 1 to 4 and the columns are indexed from 1-6 for the 4x6 plate. The desired degree of the polynomial can be specified as `spatDegree` and the names of the variables in Z that store the row and column indices can be passed in as `rowVar` and `colVar`, respectively. By default, `spatDegree=0` (indicating that no spatial effects should be included), `rowVar=:row`, and `colVar=:column`. In the example below, we have instructed `read_plate` to include spatial effects up to a second-degree polynomial. 

```
dat3 = read_plate("./data/X_df.csv", "./data/Y.csv", "./data/Zspat_df.csv", 
                  XCVar=:cond_conc, XCType="sum", 
                  ZCVar=:mut, ZCType="sum", 
                  spatDegree=2, rowVar=:row, colVar=:col)
```

Instead of reading in X, Y, and Z directly from saved files, it is possible to run `read_plate` for DataFrames that have already been loaded into the current session. The arguments to convert condition and mutant variables to contrasts, encode dosage slopes, and incorporate spatial effects remain the same. This approach may be useful if the user wants to manually pre-process the data or run `contr` on additional categorical covariates in `X` and `Z`. 

```
dat4 = read_plate(X_df, DataFrame(Y), Z_df, 
                  XCVar=:cond_conc, XCType="sum", 
                  ZCVar=:mut, ZCType="sum")
```

Finally, `read_plate` also provides the option to standardize Y by subtracting the row medians and dividing by the column IQRs if the user sets `isYstd=true`, as seen below. By default, `isYstd=false` and Y will not be standardized. 

```
dat5 = read_plate(X_df, DataFrame(Y), Z_df, 
                  XCVar=:cond_conc, XCType="sum", 
                  ZCVar=:mut, ZCType="sum", isYstd=true)
```

The `mlm` function from `matrixLM` computes least-squares coefficient estimates for matrix linear models and returns an object of type `Mlm`. `Mlm` objects include variables for the coefficient estimates (`B`), the coefficient variance estimates (`varB`), and the estimated variance of the errors (`sigma`). 

Coding the conditions in X and the mutants in Z as sum contrasts is convenient because interpretation of the main effects and interactions will be with respect to average colony sizes. However, to avoid over-parameterization, the last condition and last mutant will be left out of the contrasts fed into the model. `mlm_backest_sum` extends `mlm` by additionally back-estimating one or both of the "left-out" sum contrasts for the conditions and mutants. By default, both of the "left-out" sum contrasts will be estimated, but this behavior can be modified by setting `isXSum` and/or `isZSum` to `false`. `mlm_backest_sum` currently does not support back-estimation of contrasts when additional covariates are included in X and/or Z. 

As with `mlm`, the `mlm_backest_sum` function estimates both row and column main effects (X and Z intercepts), but this behavior can be suppressed by setting `isXIntercept=false` and/or `isZntercept=false`. Column weights for `Y` and the target type for variance shrinkage <sup>[2](#myfootnote2)</sup> can be optionally supplied to `weights` and `targetType`, respectively. 

```
est = mlm_backest_sum(dat5)
```

Note that the resulting matrix of coefficient estimates has 21 rows (20 experimental conditions plus the row intercept/main effect) and 9 columns (8 mutants plus the column intercept/main effect). 

```
size(coef(est))
```

The `coef` access function has been extended to suppress returning the "left-out" estimates if `isXSum` and/or `isZSum` is set to `true`. The `matrixLM` functions that compute predicted values (`predict`) and residuals (`resid`) have likewise been extended to include options for suppressing the "left-out" estimates. 

```
size(coef(est, isXSum=true, isZSum=true))
```

The matrix of t-statistics (defined as `est.B ./ sqrt.(est.varB)`) corresponding to this `Mlm` includes all 20x8 interactions. By default, `t_stat` does not return the corresponding t-statistics for any main effects that were estimated by `mlm_backest_sum`, but they will be returned if `isMainEff=true`. 

```
size(t_stat(est))
```

Analogous to `matrixLM`'s `mlm_perms` function, `mlm_backest_sum_perms` computes permutation p-values that include the "left-out" estimates and will run the permutations in parallel when possible. The illustrative example below only runs 5 permutations, but a different number can be specified as the second argument. By default, the function used to permute `Y` is `shuffle_rows`, which shuffles the rows for `Y`. Alternative functions for permuting `Y`, such as `shuffle_cols`, can be passed into the argument `permFun`. Keyword arguments to be passed into `mlm_backest_sum` or `matrixLM`'s `tstat` function can be specified by the user. 

```
nPerms = 5
tStats, pVals = mlm_backest_sum_perms(dat5, nPerms)
```

The `GeneticScreen` package also provides an implementation of Collins et al. (2006) <sup>[1](#myfootnote1)</sup>'s S scores. To run the `S_score` function, one must construct a RawData object where X and Z encode the experimental conditions and mutants as treatment contrasts but no intercept. This can be done by running `read_plate` with the arguments `XCType="noint"` and `ZCType="noint"`. No other covariates should be included in X and Z. 

```
dat6 = read_plate(X_df, DataFrame(Y), Z_df, 
                  XCVar=:cond_conc, XCType="noint", 
                  ZCVar=:mut, ZCType="noint", isYstd=true)
```

By default, the `S_score` function performs the variance flooring and adjustments described by Collins et al. (2006) <sup>[1](#myfootnote1)</sup> (`isVarFloor=true`). 

```
S = S_score(dat6)
```

Analogous to `matrixLM`'s `mlm_perms` function, `S_score_perms` computes permutation p-values for the S scores. 

```
S, SPvals = S_score_perms(dat6, nPerms)
```

Additional details can be found in the documentation for specific functions. 

---

<a name="myfootnote1">1</a>. Collins, S. R., Schuldiner, M., Krogan, N. J., & Weissman, J. S. (2006). A strategy for extracting and analyzing large-scale quantitative epistatic interaction data. Genome biology, 7(7), R63.

<a name="myfootnote2">2</a>. Ledoit, O., & Wolf, M. (2003). Improved estimation of the covariance matrix of stock returns with an application to portfolio selection. Journal of empirical finance, 10(5), 603-621. 