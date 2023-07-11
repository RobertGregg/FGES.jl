####################################################################
# Data structures for algorithm
####################################################################

Base.@kwdef mutable struct CurrentState
    data::Matrix{Float64}
    featureNames::Vector{String}
    numFeatures::Int
    numObservations::Int
    penalty::Float64
    
    x::Int = 1
    y::Int = 2
    bestScore::Float64 = 0.0
    bestSubset::Vector{Int} = zeros(20)
    scoreBoard::LRU{Tuple{Vector{Int}, Int}, Float64} = LRU{Tuple{Vector{Int}, Int}, Float64}(maxsize=5_000_000)
    stage::String = "Forward Search"
    verbose::Bool = true
end

function Base.show(io::IO, state::CurrentState)

    printstyled(io,"Current State\n", bold=true, color=:blue)
    println(io, "Stage: $(state.stage)")
    println(io, "Score: $(state.bestScore)")
    println(io, "Edge: $(state.x)→$(state.y)")
    println(io, "Cache: $(round(100state.scoreBoard.currentsize / state.scoreBoard.maxsize, digits=3))%")
    println(io,"----------------------------------")
end 



####################################################################
# Main entry point for the algorithm
####################################################################

"""
    fges(data, featureNames; penalty = 1.0, verbose=false, returnState=false)
Compute a causal graph for the given observed data.
"""
function fges(data, featureNames; penalty = 1.0, verbose=false, returnState=false)

    #Get the dimensions of the input data
    numObservations, numFeatures = size(data)

    #Ensure the data is in the right orientation
    if length(featureNames) ≠ numFeatures
        error("Number of features in data do not match number of feature names. Check if data has rows as observations and columns as features.")
    end

    #Append column of ones to data for intercept in linear regressions performed in scoring function
    data = [data ones(numObservations)]

    #Create an empty graph with one node for each feature
    g = SimpleDiGraph(numFeatures)

    #Create a current state for the algorithm
    state = CurrentState(
        data = data,
        featureNames = featureNames,
        numObservations = numObservations,
        numFeatures = numFeatures,
        penalty = penalty,
        verbose = verbose)

    #Forward Search
    forwardSearch!(g,state)

    #Backward Search
    state.stage = "Backward Search"
    backwardSearch!(g,state)

    if returnState
        return g, state
    end

    return g
end

#Pass only a data matrix
fges(data; args...) = fges(data, string.(1:size(data,2)); args...)

#Pass a DataFrame
fges(df::DataFrame; args...) = fges(Matrix(df), names(df); args...)


####################################################################
# Forward search
####################################################################

"""
    Insert!(g, state::CurrentState)
Modify the graph `g` by directing the edge `state.x`→`state.y`. Additionally, orient all neighbors of `y` not connected to `x` toward `y`.
"""
function Insert!(g, state::CurrentState) 

    #Add a directed edge x→y
    add_edge!(g, state.x, state.y)
    
    #Orient all edges incident into child node
    for t in state.bestSubset
        orientedge!(g, t, state.y) #t→y
    end

    return nothing
end


"""
    forwardSearch!(g, state::CurrentState)

Search equivance class space and continually add edges to `g` until the score stops increasing
"""
function forwardSearch!(g, state::CurrentState)

    #Loop until state.bestScore is negative
    while true

        #Get the next best step
        nextInsertEquivClass!(state, g)

        if state.verbose
            show(state)
        end

        if state.bestScore ≤ 0.0
            break
        end

        #Update the graph
        Insert!(g, state) 

        #Covert the PDAG to a complete PDAG
        graphVStructure!(g)

        #Apply the 4 Meek rules to orient some edges in the graph
        meekRules!(g)

        #Reset the score
        state.bestScore = 0.0
    end

    return nothing
end


"""
    nextInsertEquivClass!(state, g)
Searches through all variables pairs and finds the highest scoring edge to insert
"""
function nextInsertEquivClass!(state, g)
    
    #TODO This @sync locks the state from being updated which is probably a major bottleneck, currently needed to avoid race conditions
    @sync begin
        #Loop through all node pairs
        for (x,y) in allpairs(vertices(g))
            Threads.@spawn begin
                #Only check non-adjacent pairs
                if !isadjacent(g,x,y)
                    findBestInsert!(state, g, x, y)
                end
            end

        end
    end

    return nothing
end




"""
    findBestInsert!(state::CurrentState, g, x, y)
Scores a given variable pair and updates the state 
"""
function findBestInsert!(state::CurrentState, g, x, y)

    #The vertices neighboring y and adjacent to x need to:
        #(1) be a clique
        #(2) block all semi-directed paths
    #These will hold true regardless of T
    NAyx = calculateNAyx(g,y,x)

    if !(isclique(g, NAyx) && isblocked(g, y, x, NAyx))
        return nothing
    end


    #If those criteria are met,
    #Loop through all subsets of T
    for T in powerset(collect(calculateTyx(g,y,x)))

        #NAyx ∪ T ∪ PaY
        NTP = setdiff(ancestors(g,y), T)

        #NAyx ∪ T ∪ PaY ∪ X
        NTPx = [NTP; x] #TODO: avoid vcat in hotloop

        newScore = cached_score(state, NTPx, y) - cached_score(state, NTP, y)

        #Check if score improved
        if newScore > state.bestScore
            state.bestSubset = T
            state.bestScore = newScore
            state.x = x
            state.y = y
        end
    end

    return nothing
end

####################################################################
# Backward search
####################################################################

"""
    Delete!(g, state::CurrentState)
Modify the graph `g` by removing the edge `state.x`→`state.y`. Additionally, orient all neighbors of `x` and `y` away from `x` and `y`.
"""
function Delete!(g, state::CurrentState)

    #remove directed and unidrected edges (x→y and x-y)
    rem_edge!(g, state.x, state.y)
    rem_edge!(g, state.y, state.x)
    
    #Orient all vertices in H toward x and y
    for h ∈ state.bestSubset
        orientedge!(g, state.x, h) #x→h
        orientedge!(g, state.y, h) #y→h
    end

    return nothing
end


"""
    backwardSearch!(g, state::CurrentState)

Search equivance class space and continually remove edges in `g` until the score stops increasing
"""
function backwardSearch!(g, state::CurrentState)

    #Loop until state.bestScore is negative
    while true

        #Get the next best step
        nextDeleteEquivClass!(state, g)

        if state.verbose
            show(state)
        end

        if state.bestScore ≤ 0.0
            break
        end

        #Update the graph
        Delete!(g, state) 

        #Covert the PDAG to a complete PDAG
        graphVStructure!(g)

        #Apply the 4 Meek rules to orient some edges in the graph
        meekRules!(g)

        #Reset the score
        state.bestScore = 0.0
    end

    return nothing
end


"""
    nextInsertEquivClass!(state, g)
Searches through all variables pairs and finds the highest scoring edge to delete
"""
function nextDeleteEquivClass!(state, g)
    
    #Loop through all node pairs
    for (x,y) in allpairs(vertices(g))

        #Only check adjacent pairs
        if isadjacent(g,x,y)
            findBestDelete!(state, g, x, y)
        end
    end

    return nothing
end


"""
    findBestDelete!(state::CurrentState, g, x, y)
Scores a given variable pair and updates the state 
"""
function findBestDelete!(state::CurrentState, g, x, y)

    #Loop through all subsets of H
    for H in powerset(collect(calculateNAyx(g,y,x)))

        if isclique(g,H)

            PAy = parents(g,y)
            PAy⁻ = Iterators.filter(!isequal(x), PAy)
            newScore = cached_score(state, H ∪ PAy⁻, y) - cached_score(state, H ∪ PAy, y)

            #Check if score improved
            if newScore > state.bestScore
                state.bestSubset = H
                state.bestScore = newScore
                state.x = x
                state.y = y
            end
        end
    end

    return nothing
end


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
    
    mse = 0.0

     #If nodeParents is empty the regression is a horizontal line at mean
    if isempty(nodeParents)

        mse = sum(abs2, y .- mean(y)) / n

     #If nodeParents is of length 1, use explict formulas for simple regression
    elseif length(nodeParents) == 1

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
        mse = sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n

     #Proceed with usual linear regression
    else

        #Get the data being regression onto y
        push!(nodeParents, state.numFeatures+1)
        X  = view(state.data, :, nodeParents)

        #Perform linear regression
        b = X \ y

        #Get the estimation
        ŷ = X*b #allocations

        #Mean squared error
        mse = sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n
    end

    return mse
end