function nextInsertEquivClass!(state, g)
    
    
    (state.bestScore, state.bestSubset, state.x, state.y) = ThreadsX.mapreduce(max, allpairs(vertices(g))) do (x,y)
        if !isadjacent(g,x,y)
            findBestInsert(state, g, x, y)
        else
            (0.0,Int[],x,y)
        end
    end


    return nothing
end

function findBestInsert(state::CurrentState, g, x, y)

    bestScore = 0.0
    bestSubset = Int[]

    #The vertices neighboring y and adjacent to x need to:
        #(1) be a clique
        #(2) block all semi-directed paths
    #These will hold true regardless of T
    NAyx = calculateNAyx(g,y,x)

    if !(isclique(g, NAyx) && isblocked(g, y, x, NAyx))
        return (bestScore, bestSubset)
    end


    #If those criteria are met,
    #Loop through all subsets of T
    for T in powerset(collect(calculateTyx(g,y,x)))

        #NAyx ∪ T ∪ PaY
        NTP = setdiff(ancestors(g,y), T)

        #NAyx ∪ T ∪ PaY ∪ X
        NTPx = [NTP; x] #TOOD union! maybe?

        newScore = score(state, NTPx, y) - score(state, NTP, y)

        #Check if score improved
        if newScore > state.bestScore
            bestScore = newScore
            bestSubset = T
        end
    end

    return (bestScore, bestSubset)
end