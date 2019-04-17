# GeneticScreen

Pre- and post-processing for analysis of genetic screens. 
S-scores implemented based on Collins et al (2006) [^fn1]. This should be 
considered an extension of the `matrixLM` package. 

## Install
```
using Pkg; Pkg.add("https://github.com/janewliang/GeneticScreen.jl")
using GeneticScreen
```

[^fn1]: Collins, S. R., Schuldiner, M., Krogan, N. J., & Weissman, J. S. 
    (2006). A strategy for extracting and analyzing large-scale quantitative 
    epistatic interaction data. Genome biology, 7(7), R63.