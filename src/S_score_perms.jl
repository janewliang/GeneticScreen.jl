"""
    S_score_perms(data, nPerms; permFun, isVarFloor)

Obtains permutation p-values for S scores as described by Collins et al. 

# Arguments

- data = RawData object, where the slot for the `X` predictor matrix (row 
  covariates) is assumed to be the over-parameterized treatment contrasts for 
  the conditions and the slot for the `Z` predictor matrix (column covariates) 
  is assumed to be the over-parameterized treatment contrasts for the clones
- nPerms = number of permutations. Defaults to `1000`.

# Keyword arguments

- permFun = function used to permute `Y`. Defaults to `shuffle_rows` 
  (shuffles rows of `Y`). 
- isVarFloor = boolean indicating whether to include variance flooring and 
  adjustments. Defaults to `true`. 

# Value

Tuple
- Sscores: 2d array of floats; S scores
- pvals: 2d array of floats; permutation p-values

# Some notes

Permutations are computed in parallel when possible. 

# Reference

Collins, S. R., Schuldiner, M., Krogan, N. J., & Weissman, J. S. (2006). A 
    strategy for extracting and analyzing large-scale quantitative epistatic 
    interaction data. Genome biology, 7(7), R63.

"""
function S_score_perms(data::RawData, nPerms::Int64=1000; 
                       permFun::Function=shuffle_rows, isVarFloor::Bool=true)
    
    # Get S scores and permutation p-values
    return perm_pvals(S_score, data, nPerms; 
                      permFun=permFun, isVarFloor=isVarFloor)
end