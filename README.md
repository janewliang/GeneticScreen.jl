# GeneticScreen

Pre- and post-processing for the analysis of high-throughput genetic screens using matrix linear models. See the associated paper, ["Matrix linear models for high-throughput chemical genetic screens"](https://www.biorxiv.org/content/10.1101/468140v1), for more details. S-scores were implemented based on Collins et al. (2006) <sup>[1](#myfootnote1)</sup> . 

`GeneticScreen` is an extension of the `matrixLM` package, which provides core functions for closed-form least squares estimates for matrix linear models. 

## Install

The `GeneticScreen` package can be installed by running: 

```
using Pkg; Pkg.add("https://github.com/janewliang/GeneticScreen.jl")
```

`GeneticScreen` was developed in [Julia v1.1](https://julialang.org/downloads/). 

## Usage

```
using GeneticScreen
```

---

<a name="myfootnote1">1</a>. Collins, S. R., Schuldiner, M., Krogan, N. J., & Weissman, J. S. (2006). A strategy for extracting and analyzing large-scale quantitative epistatic interaction data. Genome biology, 7(7), R63.