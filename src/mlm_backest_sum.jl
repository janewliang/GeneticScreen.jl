"""
    backest_sum!(MLM; isXSum, isZSum)

Back-estimate left-out sum contrast interactions for an MLM object. 

# Arguments

- MLM = Mlm object

# Keyword arguments

- isXSum = boolean flag indicating whether a left-out `X` sum contrast 
  interaction needs to be back-estimated. Defaults to `true`. 
- isZSum = boolean flag indicating whether a left-out `Z` sum contrast 
  interaction needs to be back-estimated. Defaults to `true`. 

# Value

None; updates coefficient estimates in Mlm object in place. 

# Some notes

This implementation assumes that the `X` predictor matrix 
(when `isXSum` is set to `true`) and the `Z` predictor matrix 
(when `isZSum` is set to `true`) consist of the sum contrasts for a single 
categorical variable (and no other variables). 

"""
function backest_sum!(MLM::Mlm; isXSum::Bool=true, isZSum::Bool=true)

    # Check to make sure X and/or Z intercepts were included for 
    # interactions that need to be back-estimated
    if (isXSum==true) && (MLM.data.predictors.isXIntercept==false)
        error("No X intercept.")
    end
    if (isZSum==true) && (MLM.data.predictors.isZIntercept==false)
        error("No Z intercept.")
    end
    
    # Storing some larger computations to be used to calculate variances
    # Calculate and store transpose(Z)*Z
    ZTZ = transpose(get_Z(MLM.data))*get_Z(MLM.data)
    varLeft = inv(transpose(get_X(MLM.data))*get_X(MLM.data))
    varRight = transpose(ZTZ\(transpose((ZTZ\(transpose(get_Z(MLM.data))*
	                                     MLM.sigma))*get_Z(MLM.data))))

    # Back-estimate left-out X sum contrast interaction
    if isXSum==true
        C = transpose(vcat(0, -ones(MLM.data.p-1)))
        coeffX = C*MLM.B
        varX = transpose(kron_diag(C*varLeft*transpose(C), varRight))
        
        if isZSum==false
            MLM.B = vcat(MLM.B, coeffX)
            MLM.varB = vcat(MLM.varB, varX)
        end
    end
    
    # Back-estimate left-out Z sum contrast interaction
    if isZSum==true
        D = vcat(0, -ones(MLM.data.q-1))
        coeffZ = MLM.B*D
        varZ = transpose(kron_diag(varLeft, transpose(D)*varRight*D))
        
        if isXSum==false
            MLM.B = hcat(MLM.B, coeffZ)
            MLM.varB = hcat(MLM.varB, varZ)
        end
    end

    # Back-estimate left-out X and Z sum contrast interaction 
    if (isXSum==true) && (isZSum==true)
        coeffXZ = C*MLM.B*D
        varXZ = transpose(kron_diag(C*varLeft*transpose(C), 
                          transpose(D)*varRight*D))

        MLM.B = hcat(vcat(MLM.B, coeffX), vcat(coeffZ, coeffXZ))
        MLM.varB = hcat(vcat(MLM.varB, varX), vcat(varZ, varXZ))
    end
end


"""
    mlm_backest_sum(data; isXIntercept, isZIntercept, weights, targetType, 
                    isXSum, isZSum)

Matrix linear model using least squares method that optionally back-estimates 
left-out sum contrast interactions. 

# Arguments

- data = RawData object

# Keyword arguments

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
- isXSum = boolean flag indicating whether a left-out `X` sum contrast 
  interaction needs to be back-estimated. Defaults to `true`. 
- isZSum = boolean flag indicating whether a left-out `Z` sum contrast 
  interaction needs to be back-estimated. Defaults to `true`. 

# Value

None; updates coefficient estimates in Mlm object in place. 

# Some notes

This implementation assumes that the `X` predictor matrix 
(when `isXSum` is set to `true`) and the `Z` predictor matrix 
(when `isZSum` is set to `true`) consisist of the sum contrasts for a single 
categorical variable. 

"""
function mlm_backest_sum(data::RawData; 
                         isXIntercept::Bool=true, isZIntercept::Bool=true, 
                         weights=nothing, targetType=nothing, 
                         isXSum::Bool=true, isZSum::Bool=true)
    
    # Run matrix linear models
    MLM = mlm(data; isXIntercept=isXIntercept, isZIntercept=isZIntercept, 
              weights=weights, targetType=targetType)

    # Back-estimate left-out sum contrast interactions
    backest_sum!(MLM; isXSum=isXSum, isZSum=isZSum)

    return MLM
end
