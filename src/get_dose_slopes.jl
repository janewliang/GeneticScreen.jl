"""
    get_dose_slopes(X, conditionVar, concentrationVar)

Converts categorical variables in a DataFrame to dosage slopes

# Arguments

- X = DataFrame with at least two categorical variables: conditions, specified 
  by `conditionVar`, and concentrations, specified by `concentrationVar`
- conditionVar = Symbol for the categorical condition variable in `X` to use 
  for dosage slopes. 
- concentrationVar = Symbol for the categorical concentration variable in `X` 
  to use for dosage slopes. 

# Value

DataFrame

# Some notes

`get_dose_slopes` assumes that within each condition, the concentration levels 
are ordered.  

"""
function get_dose_slopes(X::DataFrames.DataFrame, 
                         conditionVar::Symbol, concentrationVar::Symbol)

    # Ensure that conditionVar is a variable in X
    if !in(conditionVar, names(X))
        error("conditionVar is not a variable in the DataFrame")
    end

    # Ensure that concentrationVar is a variable in X
    if !in(concentrationVar, names(X))
        error("concentrationVar is not a variable in the DataFrame")
    end

    # Pull out every condition 
    allConds = unique(X[:,conditionVar])
    # For each condition, pull out all of the concentrations
    concs = [X[:,concentrationVar][cond .== X[:,conditionVar]] 
             for cond in allConds]
    # Get the number of unique concentration levels found for each condition
    numLevels = [length(unique(conc)) for conc in concs]
    
    # Initialize array of slopes
    XDos = zeros(size(X, 1), length(allConds))
    
    # Iterate through conditions and assign slopes for concentrations
    for i in 1:length(allConds)
        for j in 1:length(unique(concs[i]))
            XDos[(X[:,conditionVar] .== allConds[i]) .& 
                 (X[:,concentrationVar] .== unique(concs[i])[j]), i] .= j
        end
    end
    
    # Convert results to a DataFrame and rename columns 
    XDosDf = convert(DataFrame, XDos)
    names!(XDosDf, [Symbol(cond) for cond in allConds])
    
    # Initialize new DataFrame
    newX = DataFrame()
    
	# Iterate through variables in X
    for var in names(X)
        if !in(var, [conditionVar, concentrationVar])
            # Add the other variables to the new DataFrame
            newX[!,var] = X[:,var]
        else
            # Add the DataFrame of dosage slopes		
            if var==conditionVar
                for cond in names(XDosDf)
                    newX[!,cond] = XDosDf[:,cond]
                end
            end
        end
    end
    
    return newX
end