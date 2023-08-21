using Combinatorics

#Helper function to loop through pairs of data columns
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))

function maxScore(data)

    bestScore = 0.0

    for (i,j) in allpairs(axes(data,2))
        
        if isodd(i+j)
            currentScore = findscore(data,i,j)

            if bestScore < currentScore
                bestScore = currentScore
            end
        end
    end

    return bestScore
end

function findscore(data, i, j)
    
    bestScoreSubset = 0.0

    #Create a somewhat small subset
    subsets = mod1(i+j,10)

    for subset in powerset(1:subsets)

        currentScoreSubset = score(data,subset,j)

        if bestScoreSubset < currentScoreSubset
            bestScoreSubset = currentScoreSubset
        end

    end

    return bestScoreSubset
end

function score(data,subset,j)

    if isempty(subset)
        return 0.0
    end
    
    X = view(data,:,subset)
    y = view(data,:,j)

    b = X \ y

    ŷ = X*b

    return sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ))
end



using Floops

function maxScore_par(data)

    @floop for (i,j) in allpairs(axes(data,2))

        currentScore = findscore(data,i,j)

        if isodd(i+j)
            @reduce() do (bestScore = 0.0; currentScore)
                if bestScore < currentScore
                    bestScore = currentScore
                end
            end
        end
    end

    return bestScore
end