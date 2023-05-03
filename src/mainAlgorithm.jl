#For running on a server, need to flush out the print calls
function message(x...)
    println(x...)
    flush(stdout)
end

####################################################################
# Data Structures
####################################################################

#Collection of variables to pass along to different functions in the FGES algorithm.  
struct ParseData{F<:AbstractFloat}
    data::Matrix{F} #orginal data
    normAugData::Matrix{F} #standardized by the mean and std of each column, appended with ones column at end
    scatterMat::Matrix{F} # (covariance matrix not scaled by the number of observations) 
    numFeatures::Int #number of columns
    numObservations::Int #number of rows
    ŷ::Vector{F}
    penalty::F

    function ParseData(data::Matrix{F}, scatterMat, penalty::F) where F<:AbstractFloat 

        #Get the dimensions of the input data
        numObservations, numFeatures = size(data)

        #Copy the data and standardize each column
        normAugData = copy(data)
        normAugData .-= mean(normAugData, dims=1) #subtract mean
        normAugData ./= std(normAugData, dims=1) #divide by standard deviation

        #Augment a column of ones on the end for linear regression
        normAugData = [normAugData ones(F, numObservations)]

        #If no scatter matrix is provided, calculate in construction 
        if isnothing(scatterMat)
            scatterMat = normAugData'normAugData
        end

        ŷ = zeros(F,numObservations)

        return new{F}(data, normAugData, scatterMat, numFeatures, numObservations, ŷ, penalty)
    end
end  

#Simple structure to hold the current edge, a subset of neighbors, and a score change
Base.@kwdef mutable struct Step{A,B}
    edge::Edge = Edge(1,2,true)
    subset::Vector{A} = A[]
    Δscore::B = zero(B)
end


####################################################################
# Base Function Overloads
####################################################################

#Print method to display the parsed data
function show(io::IO, d::ParseData)
    print(io, "$(d.numObservations)×$(d.numFeatures) $(typeof(d))")
end

#Print method to display the Step
function show(io::IO, newStep::Step{A,B}) where {A,B}
    print(io, "Edge: $(newStep.edge), Subset: $(newStep.subset), Δscore: $(newStep.Δscore)")
end

#This could be dangerous. The @memoize macro has to check if ParseData is the same argument. We only ever define one immutable ParseData, so it will never change. Any equality check should be true.
==(a::T, b::T) where T <: ParseData = true


####################################################################
# Main Entry point for the Algorithm
####################################################################

"""
    fges(data; debug=false)
Compute a causal graph for the given observed data.
"""
function fges(data; penalty = 1.0, scatterMat=nothing, debug=false)

    #Parse the inputted data
    dataParsed = ParseData(data, scatterMat, penalty)

    #Create an empty graph with one node for each feature
    g = PDAG(dataParsed.numFeatures)

    debug && message("Start forward search")
    #Perform the forward search 
    Search!(g, dataParsed, Insert!, debug)

    debug && message("Start backward search")
    #Perform the backward search 
    Search!(g, dataParsed, Delete!, debug)
    
    #Return the graph
    return g
end

####################################################################
# Insert and Delete Operators
####################################################################

function Insert!(g::PDAG, newStep::Step) 
    edge = newStep.edge
    T = newStep.subset

    #Add a directed edge x→y
    addedge!(g,edge)
    
    #Orient all edges incident into child node
    y = edge.child
    for t ∈ T
        orientedge!(g,t,y) #t→y
    end

    return nothing
end


function Delete!(g::PDAG, newStep::Step)
    edge = newStep.edge
    H = newStep.subset

    #remove directed and unidrected edges (x→y and x-y)
    remedge!(g, edge)
    
    #Orient all vertices in H toward x and y
    x,y = edge.parent, edge.child
    for h ∈ H
        orientedge!(g,x,h) #x→h
        orientedge!(g,y,h) #y→h
    end

    return nothing
end

####################################################################
# Forward and Backward search
####################################################################

function Search!(g, dataParsed::ParseData{A}, operator, debug) where A

    #Create a container to hold information about the next step
    newStep = Step{Int,A}()

    #Continually add/remove edges to the graph until the score stops increasing
    while true
        
        #Get the new best step (depends on if we're inserting or deleting)
        if operator == Insert!
            findNextEquivClass!(newStep, dataParsed, g, findBestInsert, debug)
        else
            findNextEquivClass!(newStep, dataParsed, g, findBestDelete, debug)
        end

        debug && statusUpdate(g)
        
        #If the score did not improve...
        if newStep.Δscore ≤ zero(A)
            #...stop searching
            break
        end
        
        #Use the insert or delete operator update the graph
        operator(g, newStep)
        
        #Covert the PDAG to a complete PDAG
        #undirect all edges unless they participate in a v-structure
        graphVStructure!(g)
        #Apply the 4 Meek rules to orient some edges in the graph
        meekRules!(g)

        #Reset the score
        newStep.Δscore = zero(A)

        
    end

    return nothing
end


function findNextEquivClass!(newStep, dataParsed, g, findBestOperation::Function, debug)

    #Parallelization with synced threads (to avoid race conditions)
    @sync begin
        #Loop through all possible node combinations
        for (x,y) in allpairs(vertices(g)) 
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
                    end
                end

            end
        end
    end
end


function findBestInsert(dataParsed::ParseData{A}, g, x, y, debug) where A
    #Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Txy: any subset of the undirected neighbors of y not connected to x
    Tyx = calcT(g, Edge(x,y,true))
    NAyx = calcNAyx(g, Edge(x,y,true))


    #Ceate two containers to hold the best found score and best subset of Tyx
    bestScore = zero(A)
    bestT = Vector{Int}()

    #Keep a list of invalid sets
    invalid = Vector{Vector{Int}}()
    
    #Loop through all possible subsets of Tyx
    for T in powerset(Tyx)
        if checkSupersets(T,invalid)
            NAyxT = NAyx ∪ T
            if isclique(g, NAyxT) && isblocked(g, y, x, NAyxT)

                #Score the valid Insert
                PAy = parents(g,y)
                PAy⁺ = PAy ∪ x
                newScore = score(dataParsed, NAyxT ∪ PAy⁺, y, debug) - score(dataParsed, NAyxT ∪ PAy, y, debug)
                
                #Save the new score if it was better than any previous
                if newScore > bestScore
                    bestT = T
                    bestScore = newScore
                end
            end
        else
            #Record that the subset T is invalid
            push!(invalid,T)
        end
    end

    return (bestT, bestScore)
end

#Check if the set T is contained in any invalid set
function checkSupersets(T,invalid)
    for i ∈ invalid
        if i ⊆ T
            return false
        end
    end
    return true
end


function findBestDelete(dataParsed::ParseData{A}, g, x, y, debug) where A
    #Calculate two (possibly empty) sets of nodes
    # NAxy: any nodes that are undirected neighbors of y and connected to x by any edge
    # Hyx: any subset of the undirected neighbors of y that are connected to x
    NAyx = calcNAyx(g, Edge(x,y,true))
    Hyx = NAyx

    #Ceate two containers to hold the best found score and best subset of Tyx
    bestScore = zero(A)
    bestH = Vector{Int}()

    #Loop through all possible subsets of Hyx
    for H in powerset(Hyx)
        #Calculate NAyx \ {H}
        NAyx_H = setdiff(NAyx,H)

        #Check if the operator is valid
        if isclique(g, NAyx_H)

            #Score the valid operator 
            PAy = parents(g,y)
            PAy⁻ = setdiff(PAy, x)
            newScore = score(dataParsed, NAyx_H ∪ PAy⁻, y, debug) - score(dataParsed, NAyx_H ∪ PAy, y, debug)
            
            if newScore > bestScore
                bestH = H
                bestScore = newScore
            end
        end
    end

    return (bestH, bestScore)
end


####################################################################
# Scoring function
####################################################################

#@memoize LRU(maxsize=10_000_000) 
@memoize LRU(maxsize=10_000_000) function score(dataParsed, nodeParents, node, debug)

    #Unpack some variables from the dataParsed structure
    n = dataParsed.numObservations #convert datatype
    data = dataParsed.data
    p = dataParsed.penalty
    #scatterMat = dataParsed.scatterMat
    ŷ = dataParsed.ŷ

    #Use views to avoid creating copies of the data in memory
        # X is the design matrix, augmented with a column of ones at the end
        # X is also been standardized so mean(columns)=0 and std(column)=1
        # y is data from the child node being tested
    @views begin
        y = data[:,node]
        X = data[:,nodeParents]
    end

    #Some small cases to check for better performance
    if isempty(nodeParents)
        ŷ .= mean(y)
    elseif length(nodeParents) == 1
        β₁, β₀ = simpleRegression(X,y)
        @. ŷ = β₁*X + β₀
    else
        #parentsAndIncept = [nodeParents; dataParsed.numFeatures+1]
        #X1 = view(data, :, parentsAndIncept)
        X1 = [X ones(n)]
        b = X1\y
        ŷ .= X1*b
    end

    # if isempty(nodeParents)
    #     parentsAndIncept = [dataParsed.numFeatures+1]
    # else
    #     parentsAndIncept = [nodeParents; dataParsed.numFeatures+1]
    # end

    # @views begin
    #     Xᵀy = scatterMat[parentsAndIncept,node]

    #     XᵀX = scatterMat[parentsAndIncept,parentsAndIncept]

    #     y = data[:,node]
    #     X = data[:,parentsAndIncept]
    # end

    # b = X \ y
    # ŷ = X*b


    #Next we want to calculate the log likelihood that these are the parents of our node
    # score = log(P(data|Model)) ≈ -BIC/2
    # because we're only comparing log likelihoods we'll ignore the 1/2 factor
    # when P(⋅) is Guassian, log(P(data|Model)) takes the form:
    # -k⋅log(n) - n⋅log(mse)
    # k is the number of free parameters and mse is mean squared error
    k = length(nodeParents)+1 #includes the intercept
    mse = sum( (yᵢ-ŷᵢ)^2 for (yᵢ,ŷᵢ) in zip(y,ŷ) ) / n

    #return the final score we want to maximize (which is technically -2BIC)
    return -p*k*log(n) - n*log(mse)
end


function simpleRegression(X,y)

    X̄ = mean(X)
    ȳ = mean(y)

    varX = 0.0
    covXy = 0.0

    for (Xᵢ, yᵢ) in zip(X,y)
        covXy += (Xᵢ-X̄)*(yᵢ-ȳ)
        varX += (Xᵢ-X̄)^2
    end

    β₁ = covXy/varX
    β₀ = ȳ - β₁*X̄
    
    return β₁, β₀
end



function statusUpdate(g::PDAG)

    numNodes = g.nv
    maxEdges = numNodes*(numNodes - 1) ÷ 2

    aveDegree = round(degreeAverage(g), sigdigits=3)

    currentTime = Dates.format(now(), dateformat"m/d/yyyy II:MM p")

    message("------------------------")
    message("Updated at $(currentTime)")

    percentEdges = round(ne(g)/maxEdges, sigdigits=3)
    message("Edges $(ne(g)) / Total $(maxEdges) ($(percentEdges)%)")
    message("Average Degree: $(aveDegree)")
end