"""
    mlm_backest_sum_perms(data; nPerms, permFun, isXIntercept, isZIntercept, 
                          weights, targetType, isXSum, isZSum, isMainEff)

Obtains permutation p-values for MLM t-statistics, including back-estimated 
left-out sum contrast interactions. 

# Arguments

- data = RawData object
- nPerms = number of permutations. Defaults to `1000`.

# Keyword arguments

- permFun = function used to permute `Y`. Defaults to `shuffle_rows` 
  (shuffles rows of `Y`). 
- isXIntercept = boolean flag indicating whether or not to include an `X` 
  intercept (row main effects). Defaults to `true`. 
- isZIntercept = boolean flag indicating whether or not to include a `Z` 
  intercept (column main effects). Defaults to `true`. 
- weights = 1d array of floats to use as column weights for `Y`, or `nothing`. 
  If the former, must be the same length as the number of columns of `Y`. 
  Defaults to `nothing`. 
- targetType = string indicating the target type toward which to shrink the 
  error variance, or `nothing`. If the former, acceptable inputs are "A", "B", 
  "C", and "D". Defaults to `nothing`.
    - "A": Target is identity matrix
    - "B": Target is diagonal matrix with constant diagonal
    - "C": Target is has same diagonal element, and same off-diagonal element
    - "D": Target is diagonal matrix with unequal entries 
- isMainEff = boolean flag indicating whether or not to include p-values for 
  the main effects
- isXSum = boolean flag indicating whether an left-out `X` sum contrast 
  interaction needs to be back-estimated
- isZSum = boolean flag indicating whether an left-out `Z` sum contrast 
  interaction needs to be back-estimated

# Value

Tuple
- tStats: 2d array of floats; t-statistics
- pvals: 2d array of floats; permutation p-values

# Some notes

Permutations are computed in parallel when possible. 

This implementation assumes that the `X` predictor matrix 
(when `isXSum` is set to `true`) and the `Z` predictor matrix 
(when `isZSum` is set to `true`) consist of the sum contrasts for a single 
categorical variable (and no other variables). 

"""
function mlm_backest_sum_perms(data::RawData, nPerms::Int64=1000; 
                               permFun::Function=shuffle_rows, 
                               isXIntercept::Bool=true, 
                               isZIntercept::Bool=true, 
                               weights=nothing, targetType=nothing, 
                               isXSum::Bool=true, isZSum::Bool=true, 
                               isMainEff::Bool=false) 
    
    # Wrapper function that performs MLM, back-estimates left-out sum sum 
    # contrast interactions, and gets t-statistics
    function mlm_backest_sum_t_stat(data::RawData)
        return t_stat(mlm_backest_sum(data; isXIntercept=isXIntercept, 
                                      isZIntercept=isZIntercept, 
                                      weights=weights, targetType=targetType, 
                                      isXSum=isXSum, isZSum=isZSum), isMainEff)
    end
    
    # Get MLM t-statistcs and permutation p-values
    return perm_pvals(mlm_backest_sum_t_stat, data, nPerms; permFun=permFun)
end