"""
    coef(MLM; isXSum, isZSum)

Extracts coefficients from Mlm object, accounting for the possibility of the 
object containing back-estimated interactions

# Arguments 

- MLM = Mlm object

# Keyword arguments

- isXSum = boolean flag indicating if a left-out `X` sum contrast 
  interaction was back-estimated AND should not be included in the extracted 
  coefficients. Defaults to `false`. 
- isZSum = boolean flag indicating if a left-out `Z` sum contrast 
  interaction was back-estimated AND should not be included in the extracted 
  coefficients. Defaults to `false`. 

# Value

2d array of floats

# Some notes

If the left-out `X` and/or `Z` sum contrast interactions were back-estimated, 
but the user wishes to include them in the extracted coefficients anyway, set 
`isXSum` and/or `isZSum` to `false`. 

"""
function coef(MLM::Mlm; isXSum::Bool=false, isZSum::Bool=false)
    
    return MLM.B[1:(end-isXSum), 1:(end-isZSum)]
end


"""
    predict(MLM, newPredictors; isXSum, isZSum)

Calculates new predictions based on Mlm object, accounting for the possibility 
of the object containing back-estimated interactions

# Arguments 

- MLM = Mlm object
- newPredictors = Predictors object. Defaults to the `data.predictors` field 
  in the Mlm object used to fit the model. 

# Keyword arguments

- isXSum = boolean flag indicating whether a left-out `X` sum contrast 
  interaction was back-estimated. Defaults to `false`. 
- isZSum = boolean flag indicating whether a left-out `Z` sum contrast 
  interaction was back-estimated. Defaults to `false`. 

# Value

Response object

"""
function predict(MLM::Mlm, newPredictors::Predictors=MLM.data.predictors; 
                 isXSum::Bool=false, isZSum::Bool=false)
    
  	# Include X and Z intercepts in new data if necessary
  	if MLM.data.predictors.isXIntercept==true && 
       newPredictors.isXIntercept==false
    	newPredictors.X = add_intercept(newPredictors.X)
    	newPredictors.isXIntercept = true
    	println("Adding X intercept to newPredictors.")
  	end
  	if MLM.data.predictors.isZIntercept==true && 
       newPredictors.isZIntercept==false
    	newPredictors.Z = add_intercept(newPredictors.Z)
    	newPredictors.isZIntercept = true
    	println("Adding Z intercept to newPredictors.")
  	end
    
  	# Remove X and Z intercepts in new data if necessary
  	if MLM.data.predictors.isXIntercept==false && 
       newPredictors.isXIntercept==true
    	newPredictors.X = remove_intercept(newPredictors.X)
    	newPredictors.isXIntercept = false
    	println("Removing X intercept from newPredictors.")
  	end
  	if MLM.data.predictors.isZIntercept==false && 
       newPredictors.isZIntercept==true
    	newPredictors.Z = remove_intercept(newPredictors.Z)
    	newPredictors.isZIntercept = false
    	println("Removing Z intercept from newPredictors.")
  	end
    
    # Calculate predictions
    return Response(calc_preds(newPredictors.X, newPredictors.Z, 
                               coef(MLM, isXSum=isXSum, isZSum=isZSum))) 
end 


"""
    fitted(MLM; isXSum, isZSum)

Calculate fitted values of an Mlm object, accounting for the possibility of 
the object containing back-estimated interactions

# Arguments 

- MLM = Mlm object

# Keyword arguments

- isXSum = boolean flag indicating whether a left-out `X` sum contrast 
  interaction was back-estimated. Defaults to `false`. 
- isZSum = boolean flag indicating whether a left-out `Z` sum contrast 
  interaction was back-estimated. Defaults to `false`. 


# Value

Response object

"""
function fitted(MLM::Mlm; isXSum::Bool=false, isZSum::Bool=false)
    
    # Call the predict function with default newPredictors
    return predict(MLM, isXSum=isXSum, isZSum=isZSum)
end


"""
    resid(MLM, newData; isXSum, isZSum)

Calculates residuals of an Mlm object, accounting for the possibility of the 
object containing back-estimated interactions

# Arguments 

- MLM = Mlm object
- newData = RawData object. Defaults to the `data` field in the Mlm object 
  used to fit the model. 

# Keyword arguments

- isXSum = boolean flag indicating whether a left-out `X` sum contrast 
  interaction was back-estimated. Defaults to `false`. 
- isZSum = boolean flag indicating whether a left-out `Z` sum contrast 
  interaction was back-estimated. Defaults to `false`. 

# Value

2d array of floats

"""
function resid(MLM::Mlm, newData::RawData=MLM.data; 
               isXSum::Bool=false, isZSum::Bool=false)
    
    # Include X and Z intercepts in new data if necessary
    if MLM.data.predictors.isXIntercept==true && 
       newData.predictors.isXIntercept==false
        newData.predictors.X = add_intercept(newData.predictors.X)
        newData.predictors.isXIntercept = true
        newData.p = newData.p + 1
        println("Adding X intercept to newData.")
    end
    if MLM.data.predictors.isZIntercept==true && 
       newData.predictors.isZIntercept==false
        newData.predictors.Z = add_intercept(newData.predictors.Z)
        newData.predictors.isZIntercept = true
        newData.q = newData.q + 1
        println("Adding Z intercept to newData.")
    end
    
    # Remove X and Z intercepts in new data if necessary
    if MLM.data.predictors.isXIntercept==false && 
       newData.predictors.isXIntercept==true
        newData.predictors.X = remove_intercept(newData.predictors.X)
        newData.predictors.isXIntercept = false
        newData.p = newData.p - 1
        println("Removing X intercept from newData.")
    end
    if MLM.data.predictors.isZIntercept==false && 
       newData.predictors.isZIntercept==true
        newData.predictors.Z = remove_intercept(newData.predictors.Z)
        newData.predictors.isZIntercept = false
        newData.q = newData.q - 1
        println("Removing Z intercept from newData.")
    end
    
	  # Calculate residuals
    return calc_resid(get_X(newData), get_Y(newData), get_Z(newData), 
                      coef(MLM, isXSum=isXSum, isZSum=isZSum)) 
end
