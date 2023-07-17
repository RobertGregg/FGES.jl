####################################################################
# Main entry point for the algorithm
####################################################################

"""
    fges(data, featureNames; penalty = 1.0, verbose=false, returnState=false)
Compute a causal graph for the given observed data.
"""
function fges(data, featureNames; penalty = 1.0, verbose=false, returnState=false)

    
    #Ensure the data is in the right orientation
    if length(featureNames) ≠ size(data,2)
        error("Number of features in data do not match number of feature names. Check if data has rows as observations and columns as features.")
    end
    
    
    #Create a current state for the algorithm
    state = CurrentState(data, featureNames, penalty, verbose)
    
    #Create an empty graph with one node for each feature
    g = SimpleDiGraph(state.numFeatures)

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


#TODO Need to overhaul this parallelization.
#=
    1) Make a serial and parallel version to experiment. 
    2) LinearSolveBuffer needs to be local on each thread
        - How do you avoid re-initialization each loop?
        - https://julialang.org/blog/2023/07/PSA-dont-use-threadid/ 
    3) R doesn't need to be copied to each thread
    4) nextInsertEquivClass! and the like should be pure functions (currently LinearSolveBuffer and state are modified) 
    5) Is the locking of LRUCache a bottleneck?

=#

"""
    nextInsertEquivClass!(state, g)
Searches through all variables pairs and finds the highest scoring edge to insert
"""
function nextInsertEquivClass!(state, g)
    
    for (x,y) in allpairs(vertices(g))
        if !isadjacent(g,x,y)
            findBestInsert(state, g, x, y)
        end
    end

    return nothing
end


"""
    findBestInsert!(state::CurrentState, g, x, y)
Scores a given variable pair and updates the state 
"""
function findBestInsert(state::CurrentState, g, x, y)

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
        NTPx = [NTP; x]

        newScore = cached_score(state, NTPx, y) - cached_score(state, NTP, y)


        #Check if score improved
        if newScore > state.bestScore
            state.bestScore = newScore
            state.bestSubset = T
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