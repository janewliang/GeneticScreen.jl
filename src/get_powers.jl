"""
    get_powers(x, y, powers)

Obtain a DataFrame of the powers and cross products of x and y 

# Arguments

- x = 1d array of floats
- y = 1d array of floats of same length as `x`
- powers = 1d array of integers indicating the powers to include. Defaults to 
  [1, 2, 3, 4] to get everything up to the quadric term

# Value

DataFrame of the powers and cross products of `x` and `y`, as specified by the 
powers argument. 

"""
function get_powers(x::Array{Float64,1}, y::Array{Float64,1}, 
                    powers::Array{Int64,1}=[1, 2, 3, 4])
    
    # Check that x and y are the same length
    @assert length(x) == length(y)

    # Initialze DataFrame to return
    out = DataFrame() 
    
    # Get powers of x and y
    for i in 1:length(powers)
        out[:, Symbol(string("X", powers[i]))] = x.^powers[i]
            out[:, Symbol(string("Y", powers[i]))] = y.^powers[i]
    end
    
    # Get cross products
    for i in 1:(length(powers)-1)
        for j in 1:(length(powers)-i)
            out[!, Symbol(string("X", powers[i], "Y", powers[j]))] = 
                out[:, Symbol(string("X", powers[i]))] .* 
                out[:, Symbol(string("Y", powers[j]))]
        end
    end
    
    return out
end
