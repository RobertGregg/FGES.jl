####################################################################
# Scoring function
####################################################################

"""
    score(state::CurrentState, nodeParents, node)

Calculates a score using the mean-squared error found by regessing `node` onto the `nodeParents`. 

More Info: we want to calculate the log likelihood that `nodeParents` are the parents of our node

    score = log(P(data|Model)) ≈ -BIC/2

because we're only comparing log likelihoods we'll ignore the 1/2 factor. When P(⋅) is Guassian, log(P(data|Model)) takes the form:
    score = -k⋅log(n) - n⋅log(mse)
k is the number of free parameters, n is the number of observations, and mse is mean squared error
"""
function score(state::CurrentState, nodeParents, node)
    
    #Number of observations
    n = state.numObservations
    #Penalty value
    p = state.penalty
    
    #Get the data for the node to be scored
    y = view(state.data, :, node)
    
    #Calculate the mean squared error
    mse = calculateMSE(state, nodeParents, y, n)
    
    #Number of free parameters (includes numFeatures+1)
    k = length(nodeParents)

    newscore = -p*k*log(n) - n*log(mse)

    #Remove numFeatures+1 if needed
    if state.numFeatures+1 ∈ nodeParents
        pop!(nodeParents)
    end

    #Return the score
    return newscore
end


function cached_score(state::CurrentState, nodeParents, node)
    get!(state.scoreBoard, (nodeParents, node)) do
        score(state, nodeParents, node)
    end
end

"""
    calculateMSE(state, nodeParents, y, n)
Calculate the mean square error dependent on the length of `nodeParents`.
"""
function calculateMSE(state, nodeParents, y, n)
    

    #If nodeParents is empty the regression is a horizontal line at mean
    if isempty(nodeParents)
        return sum(abs2, y .- mean(y)) / n
    end


    #If nodeParents is of length 1, use explict formulas for simple regression
    if length(nodeParents) == 1

        #ŷ = β₁X + β₀
        #β₁ = cov(X,y)/var(X)
        #β₀ = mean(y) - β₁ * mean(X)

        #Get the X vector
        X  = view(state.data, :, nodeParents)

        X̄ = mean(X)
        ȳ = mean(y)

        varX = 0.0
        covXy = 0.0

        for (Xᵢ, yᵢ) in zip(X,y)
            covXy += (Xᵢ-X̄)*(yᵢ-ȳ)
            varX += (Xᵢ-X̄)^2
        end

        β₁ = covXy/varX
        β₀ = ȳ - β₁*X̄
    
        ŷ = @. β₁ * X + β₀

        #Mean squared error
        return sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n
    end

    #Proceed with usual linear regression

    #Get the data being regression onto y
    push!(nodeParents, state.numFeatures+1)

    @views begin
        X  = state.data[:, nodeParents]
        b = state.buffer.b[1:length(nodeParents)]
        ŷ = state.buffer.ŷ
    end

    #Perform linear regression
    #b = X \ y
    qless!(b, X, y, state.buffer.R, state.buffer.Rsub, nodeParents)

    #Get the estimation
    #ŷ = X*b
    mul!(ŷ, X, b)

    #Mean squared error
    return sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n
end


####################################################################
# Solve X\y using "Q-less" QR decomposition
####################################################################

#=
    X*b = y         multiply by X'
    X'X*b = X'y     Subsitute X=QR noting rules of transposes: X'=R'Q'

    R'Q'QR*b = X'y  Q is orthogonal meaning Q'Q = I (identity matrix)
    R'R*b = X'y

 This equation is useful because you only have to calculate R once. Then if you want to solve x\y where x is a subset of X, you can efficiently subset R and recalculate b.
=#

#Loop through lower triangular matrix
#the 1st axis is reversed because we're moving non-zeros *up* above diagonal
lowTriIter(A::AbstractMatrix) = ((i,j) for i in reverse(axes(A,1)), j in axes(A,2) if i>j)

function qless!(b,X,y,R,Rsub,subset)

    numCol = length(subset)

    #Rsub = copy(R[:,subset]) #dont want to destroy R
    #copyto!(Rsub, R[:,subset])
    for (k,j) in enumerate(subset)
        @views Rsub[:,k] = R[:,j]
    end

    #Use a givens rotation to zero out values below diagonal
    @inbounds for (i,j) in lowTriIter(Rsub)        
        if Rsub[i,j] ≠ 0.0
            G,r = givens(Rsub,i-1,i,j)
            Rsub[i-1,j] = r
            Rsub[i,j] = 0.0

            for k in j+1:numCol
                (r1, r2) = Rsub[i-1,k], Rsub[i,k]
                Rsub[i-1,k] = G.c*r1 + G.s*r2
                Rsub[i,k] = -G.s*r1 + G.c*r2
            end
        end
    end

    Rv = UpperTriangular(view(Rsub,1:numCol,1:numCol))

    #Solve R'R*b = X'y
    mul!(b,X',y)
    ldiv!(Rv',b)
    ldiv!(Rv,b)

    return nothing
end