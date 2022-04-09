#Using FLoops to parallelize findNextEquivClass!
using FLoops

#Generate an iterator equivalent to combinations(x,2) but with base packages
#The combinations are not in order, but that doesn't matter in this case
#benefit: packages have more methods defined for the base iterators
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))

function findNextEquivClass!(newStep, dataParsed, g, findBestOperation::Function, debug)

    #Parallelization with synced threads (to avoid race conditions)
    @sync begin
        #Loop through all possible node combinations
        for (x,y) in combinations(vertices(g),2) 
            #Spawn a new thread
            Threads.@spawn begin
                #Check if there is an edge between two vertices
                hasEdge = isadjacent(g,x,y)

                #Skip over adjacent edges if we're trying to insert and vice versa
                if findBestOperation == findBestInsert ? !hasEdge : hasEdge
                    #For this pair of nodes, the following function will:
                        #(1) check if a valid operator exists
                        #(2) score all valid operators and return the one with the highest score
                        #findBestOperation is either "findBestInsert" or "findBestDelete"
                    (bestSubset, bestScore) = findBestOperation(dataParsed, g, x, y, debug)

                    #If the best valid operator was better than any previously found...
                    if bestScore > newStep.Δscore 
                        #...update newStep
                        newStep.edge = Edge(x,y,true)
                        newStep.subset = bestSubset
                        newStep.Δscore = bestScore

                        debug && message(newStep)
                    end
                end
            end
        end
    end
end

ThreadsX.maximum((f(x,y),x,y) for (x,y) in allpairs(v))

maximum(allpairs(v)) do (x,y)
    (f(x,y),x,y)
end



function testPar!(newStep, dataParsed, g, findBestOperation::Function, debug)

    for (x,y) in allpairs(vertices(g))
        
    end
    
end