"""
    read_plate(X, Y, Z; 
               XCVar, ZCVar, XCType, ZCType, 
               isXDos, XconditionVar, XconcentrationVar, isYstd, 
               spatDegree, rowVar, colVar)

Convert the `X`, `Y`, and `Z` arrays corresponding to a plate to a RawData 
object. Optionally: 

- Specify contrasts for `X` and/or `Z`, including the option for encoding `X` 
  as dosage slopes
- Standardize the rows of the response matrix `Y` by median and IQR
- Calculate spatial effects based on row/column variables in the `Z` matrix, 
  up to a specified polynomial degree

# Arguments

- X = DataFrame with the `X` predictor matrix (row covariates)
- Y = DataFrame with the multivariate `Y` response matrix
- Z = DataFrame with the `Z` predictor matrix (column covariates)

# Keyword arguments

- XCVar = Symbol for the categorical variable in `X` to be converted into 
  dummy indicators for the conditions. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- ZCVar = Symbol for the categorical variable in `Z` to be converted into 
  dummy indicators for the mutants. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- XCType = character string indicating the type of contrast to use for 
  `XCVar`. Defaults to treatment contrasts ("treat"). Other options include 
  "sum" for sum contrasts, "noint" for treatment contrasts with no intercept, 
  and "sumnoint" for sum contrasts with no intercept. 
- ZCType = character string indicating the type of contrast to use for 
  `ZCVar`. Defaults to treatment contrasts ("treat"). Other options include 
  "sum" for sum contrasts, "noint" for treatment contrasts with no intercept, 
  and "sumnoint" for sum contrasts with no intercept. 
- isXDos = boolean flag indicating whether the `X` matrix should be converted 
  to dosage slopes. Defaults to `false`. 
- XConditionVar = Symbol for the categorical condition variable in `X` to use 
  for dosage slopes. To be passsed into `get_dose_slopes`. Defaults to empty 
  Symbol, `Symbol()`. 
- XConcentrationVar = Symbol for the categorical concentration variable in `X` 
  to use for dosage slopes. To be passsed into `get_dose_slopes`. Defaults to 
  empty Symbol, `Symbol()`. 
- isYstd = boolean flag indicating whether to standardize rows of `Y` by 
  median and IQR. Defaults to `false`. 
- spatDegree = polynomial degree for spatial effects. Defaults to 0, meaning 
  that no spatial effects will be calculated. If `spatDegree` is greater than 
  0, spatial effects will be calculated based on `rowVar` and `colVar`. 
  `rowVar` and `colVar` will then be dropped from `Z`. 
- rowVar = name of variable in Z that stores the plate rows, to be used to 
  calculate spatial effects when `spatDegree` is greater than 0. Defaults to 
  `:row`. 
- colVar = name of variable in Z that stores the plate columns, to be used to 
  calculate spatial effects when `spatDegree` is greater than 0. Defaults to 
  `:column`. 

# Value

RawData object

# Some notes

This is the base `read_plate` function. 

"""
function read_plate(X::DataFrames.DataFrame, Y::DataFrames.DataFrame, 
                    Z::DataFrames.DataFrame; 
                    XCVar::Symbol=Symbol(), ZCVar::Symbol=Symbol(),
                    XCType::String="treat", ZCType::String="treat", 
                    isXDos::Bool=false, 
                    XConditionVar::Symbol=Symbol(), 
                    XConcentrationVar::Symbol=Symbol(), 
                    isYstd::Bool=false, spatDegree::Int64=0, 
                    rowVar::Symbol=:row, colVar::Symbol=:column)
    
    # Convert Y to an array
    Y = convert(Array{Float64, 2}, Y)

    # Standardize rows of Y by median and IQR
    if isYstd == true 
        for i in 1:size(Y,1)
            Y[i,:] = (Y[i,:] .- median(vec(Y[i,:]))) ./ iqr(vec(Y[i,:]))
        end
    end
    
    # If spatDegree is specified as a positive integer, add spatial effects 
    # and delete the row/column variables from Z
    if spatDegree > 0
        centerRow = mean([maximum(Z[!, rowVar]), minimum(Z[!, rowVar])])
        centerCol = mean([maximum(Z[!, colVar]), minimum(Z[!, colVar])])
        Z = hcat(Z, get_powers(Z[!, rowVar].-centerRow, 
                               Z[!, colVar].-centerCol, 
                               collect(1:spatDegree)))
        select!(Z, Not(rowVar))
        select!(Z, Not(colVar))
    end
    
    # Create contrasts for the X and Z covariates and convert them to arrays
    if isXDos == true
        XContr = convert(Array{Float64, 2}, get_dose_slopes(X, XConditionVar, 
                                                            XConcentrationVar))
    else
        XContr = convert(Array{Float64, 2}, contr(X, [XCVar], [XCType]))
    end
    ZContr = convert(Array{Float64, 2}, contr(Z, [ZCVar], [ZCType]))
    
    # Convert arrays to RawData object
    return RawData(Response(Y), Predictors(XContr, ZContr))
end


"""
    read_plate(predictorXPath, responseYPath, predictorZPath; 
               XCVar, ZCVar, XCType, ZCType, 
               isXDos, XconditionVar, XconcentrationVar, isYstd, 
               spatDegree, rowVar, colVar, funArgs...)

Read in the `X`, `Y`, and `Z` arrays corresponding to a plate and convert 
them into a RawData object. Optionally: 

- Specify contrasts for `X` and/or `Z`, including the option for encoding `X` 
  as dosage slopes
- Standardize the rows of the response matrix `Y` by median and IQR
- Calculate spatial effects based on row/column variables in the `Z` matrix, 
  up to a specified polynomial degree

# Arguments

- predictorXPath = string with path to flat file storing the `X` predictor 
  matrix (row covariates)
- responseYPath = string with path to flat file storing the multivariate `Y` 
  response matrix
- predictorZPath = string with path to flat file storing the `Z` predictor 
  matrix (column covariates)

# Keyword arguments

- XCVar = Symbol for the categorical variable in `X` to be converted into 
  dummy indicators for the conditions. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- ZCVar = Symbol for the categorical variable in `Z` to be converted into 
  dummy indicators for the mutants. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- XCType = character string indicating the type of contrast to use for 
  `XCVar`. Defaults to treatment contrasts ("treat"). Other options include 
  "sum" for sum contrasts, "noint" for treatment contrasts with no intercept, 
  and "sumnoint" for sum contrasts with no intercept. 
- ZCType = character string indicating the type of contrast to use for 
  `ZCVar`. Defaults to treatment contrasts ("treat"). OOther options include 
  "sum" for sum contrasts, "noint" for treatment contrasts with no intercept, 
  and "sumnoint" for sum contrasts with no intercept. 
- isXDos = boolean flag indicating whether the `X` matrix should be converted 
  to dosage slopes. Defaults to `false`. 
- XConditionVar = Symbol for the categorical condition variable in `X` to use 
  for dosage slopes. To be passsed into `get_dose_slopes`. Defaults to empty 
  Symbol, `Symbol()`. 
- XConcentrationVar = Symbol for the categorical concentration variable in `X` 
  to use for dosage slopes. To be passsed into `get_dose_slopes`. Defaults to 
  empty Symbol, `Symbol()`. 
- isYstd = boolean flag indicating whether to standardize rows of `Y` by 
  median and IQR. Defaults to `false`. 
- spatDegree = polynomial degree for spatial effects. Defaults to 0, meaning 
  that no spatial effects will be calculated. If `spatDegree` is greater than 
  0, spatial effects will be calculated based on `rowVar` and `colVar`. 
  `rowVar` and `colVar` will then be dropped from `Z`. 
- rowVar = name of variable in Z that stores the plate rows, to be used to 
  calculate spatial effects when `spatDegree` is greater than 0. Defaults to 
  `:row`. 
- colVar = name of variable in Z that stores the plate columns, to be used to 
  calculate spatial effects when `spatDegree` is greater than 0. Defaults to 
  `:column`. 
- funArgs = variable keyword arguments to be passed into `CSV.read`

# Value

RawData object

# Some notes

Files are read in using `CSV.read`. By default, the first row is assumed to be 
a header and the delimiter is taken to be ','.

"""
function read_plate(predictorXPath::String, responseYPath::String, 
                    predictorZPath::String; 
                    XCVar::Symbol=Symbol(), ZCVar::Symbol=Symbol(),
                    XCType::String="treat", ZCType::String="treat", 
                    isXDos::Bool=false, 
                    XConditionVar::Symbol=Symbol(), 
                    XConcentrationVar::Symbol=Symbol(), 
                    isYstd::Bool=false, spatDegree::Int64=0, 
                    rowVar::Symbol=:row, colVar::Symbol=:column, funArgs...)
    
    # Read in the response and predictor arrays
    X = CSV.read(predictorXPath; funArgs...)
    Y = CSV.read(responseYPath; funArgs...)
    Z = CSV.read(predictorZPath; funArgs...)
    
    # Pass the arrays to the base read_plate function
    return read_plate(X, Y, Z; 
                      XCVar=XCVar, ZCVar=ZCVar, XCType=XCType, ZCType=ZCType, 
                      isXDos=isXDos, XConditionVar=XConditionVar, 
                      XConcentrationVar=XConcentrationVar, isYstd=isYstd, 
                      spatDegree=spatDegree, rowVar=rowVar, colVar=colVar)
end