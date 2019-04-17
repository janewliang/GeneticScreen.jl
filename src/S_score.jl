"""
    S_score(data; isVarFloor)

Function to get S scores as described by Collins et al. Optionally 
includes variance flooring and adjustments. 

# Arguments

- data = RawData object, where the slot for the `X` predictor matrix 
  (row covariates) is assumed to be the over-parameterized treatment  
  contrasts for the conditions and the slot for the `Z` predictor matrix 
  (column covariates) is assumed to be the over-parameterized treatment 
  contrasts for the clones
- isVarFloor = boolean indicating whether to include variance flooring and 
  adjustments. Defaults to `true`. 

# Value

2d array of floats

# Reference

Collins, S. R., Schuldiner, M., Krogan, N. J., & Weissman, J. S. (2006). A 
    strategy for extracting and analyzing large-scale quantitative epistatic 
    interaction data. Genome biology, 7(7), R63.

"""
function S_score(data::RawData; isVarFloor::Bool=true)
    
    # Pull out arrays from RawData object 
    # Response matrix with colony sizes/opacities
    Y = get_Y(data) 
    # Over-parameterized treatment contrasts for the conditions
    allConds = get_X(data) 
    # Over-parameterized treatment contrasts for the clones
    allClones = get_Z(data) 
    
    # Allocate array to store S scores. 
    S = Array{Float64}(undef, size(allConds,2), size(allClones,2))
    
    # Allocate array to store mean/median colony size over clones
    muCont = Array{Float64}(undef, size(allClones,2))
    
    # If isVarFloor is true: get median colony size over clones
    if (isVarFloor == true) 
        for j=1:size(allClones,2) # Iterate through all clones
            thisClone = Y[:, allClones[:,j].==1]

            # median(colony size for that clone)
            muCont[j] = Statistics.median(thisClone) 
        end
    # Otherwise, get mean, variance, and number of measurements of colony size 
    # over clones
    else 
        # Allocate arrays for storing colony size variance and number of 
        # measurements of colony sizes over clones
        varCont = Array{Float64}(undef, size(allClones,2))
        nCont = Array{Float64}(undef, size(allClones,2))
        for j=1:size(allClones,2) # Iterate through all clones
            thisClone = Y[:, allClones[:,j].==1]

            # mean(colony size for that clone)
            muCont[j] = Statistics.mean(thisClone) 
            # var(colony size for that clone)
            varCont[j] = Statistics.var(thisClone) 
            # num of measurements of colony sizes for that clone
            nCont[j] = length(thisClone) 
        end
    end 
    
    # Allocate arrays to store mean, variance, and number of measurements of 
    # colony sizes over clones and conditions
    muExp = Array{Float64}(undef, size(allConds,2), size(allClones,2)) 
    varExp = Array{Float64}(undef, size(allConds,2), size(allClones,2))
    nExp = Array{Float64}(undef, size(allConds,2), size(allClones,2)) 
    
    # Iterate through all conditions and clones
    for i=1:size(allConds,2), j=1:size(allClones,2) 
        
        thisCloneCond = Y[allConds[:,i].==1, allClones[:,j].==1]

        # mean(colony sizes for that clone and condition)
        muExp[i,j] = Statistics.mean(thisCloneCond) 
        # var(colony sizes for that clone and condition)
        varExp[i,j] = Statistics.var(thisCloneCond) 
        # num of measurements of colony sizes for that clone and condition
        nExp[i,j] = length(thisCloneCond) 
    end
    
    # If isVarFloor is true: control the lower bounds for variance and get 
    # the S score
    if (isVarFloor == true) 
        # Control variance with lower bound. 
        varContLB = (muCont * Statistics.median(sqrt.(varExp)./muExp)).^2
        varCont = max.(transpose(Statistics.median(varExp, dims=1)), varContLB)
        
        # Control sample size. 
        nCont = Statistics.median(nExp)
        
        # Experimental variance with lower bound. 
        sdExpLoess = loess(vec(muExp), vec(sqrt.(varExp)))
        varExp = reshape(max.(vec(varExp), 
                              Loess.predict(sdExpLoess, vec(muExp)).^2), 
                         size(muExp))
        
        # Calculate S score. 
        sVar = (varExp.*(nExp.-1) .+ transpose(varCont.*(nCont.-1))) ./ 
               (nExp .+ nCont .- 2)
        S = (muExp.-transpose(muCont)) ./ sqrt.(sVar./nExp.+sVar./nCont)
        
    else # Otherwise, just get the S score
        # Calculate S score. 
        S = (muExp.-transpose(muCont)) ./ 
            sqrt.(varExp./nExp .+ transpose(varCont./nCont))
    end 
    
    return S
end