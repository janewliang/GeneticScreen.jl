module GeneticScreen


using matrixLM
import matrixLM.calc_preds, matrixLM.calc_preds!, 
    matrixLM.calc_resid, matrixLM.calc_resid!

using DataFrames
using Loess
using Statistics
import StatsBase.iqr
using CSV


export Response, Predictors, RawData, get_X, get_Z, get_Y, 
    add_intercept, remove_intercept, shuffle_rows, shuffle_cols, contr, 
    Mlm, mlm, t_stat, perm_pvals, mlm_perms, 
    get_powers, get_dose_slopes, read_plate, 
    backest_sum!, mlm_backest_sum, mlm_backest_sum_perms, 
    coef, predict, fitted, resid, 
    S_score, S_score_perms


# Get powers and cross products (for spatial effects)
include("get_powers.jl")

# Dosage slopes
include("get_dose_slopes.jl") 

# Read in plate data and convert to RawData object
include("read_plate.jl")

# MLM that back-estimates sum contrasts
include("mlm_backest_sum.jl")
# Permutations
include("mlm_backest_sum_perms.jl")
# Predictions and residuals
include("predict.jl")

# Collins S score and permutations
include("S_score.jl")
include("S_score_perms.jl")


end 
