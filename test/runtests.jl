using Revise
using FGES
using Combinatorics
using Memoization, LRUCache #Caching scores
using CSV, DataFrames
using Test


#Generate a dataset
numFeatures = 50
numObservations = 10000
data = zeros(numObservations, numFeatures)

for i in 1:numFeatures
    if i ≤ 2
        data[:,i] = randn(numObservations)
    else
        data[:,i] = sum(rand()*data[:,i-j] for j∈1:2) + randn(numObservations)
    end
end

g = fges(data)



#Save the dataset to run for java
df = DataFrame(data, :auto)
CSV.write("test/testDatasets/features$(numFeatures)_observations$(numObservations).txt", df ,delim=" ")


# #Run the julia FGES and time it
Memoization.empty_all_caches!()
g = @time fges(data)

@profview fges(data)

#Read and parse the java results
gJava = PDAG(numFeatures)
javaFilepath = "test/testDatasets/features$(numFeatures)_observations$(numObservations)_java.txt"

javaResult = filter(line -> occursin(r"x\d+ \D{3} x\d+",line),readlines(open(javaFilepath)))
javaResult = split.(javaResult," ")

javaEdges = [parse.(Int, filter.(isdigit, line[[2,4]])) for line in javaResult]

for edge in javaEdges
    addedge!(gJava,edge...,directed=false) #not worried about orientation for now
end


#Generate the true graph
gTrue = PDAG(numFeatures)

for i=3:numFeatures
    addedge!(gTrue,i-1, i)
    addedge!(gTrue,i-2, i)
end


#Calculate precision and recall
function contingencyTable(g,gTrue)

    continTable = zeros(Int,2,2)

    for (x,y) in combinations(vertices(gTrue),2)
            trueEdge = isadjacent(gTrue,x,y)
            estEdge = isadjacent(g,x,y)

            if trueEdge & estEdge #true positive
                continTable[1,1] += 1
            elseif !trueEdge & !estEdge #true negative
                continTable[2,2] += 1
            elseif !trueEdge & estEdge #false positive
                continTable[2,1] += 1
            else #false negative
                continTable[1,2] += 1
            end
    end

    return continTable
end
 
modelPrecision(continTable) = continTable[1,1] / sum(continTable[:,1])
modelRecall(continTable) = continTable[1,1] / sum(continTable[1,:])


continTable = contingencyTable(g,gTrue)
modelPrecision(continTable)
modelRecall(continTable)

continTable = contingencyTable(gJava,gTrue)
modelPrecision(continTable)
modelRecall(continTable)



g = PDAG(8)
addedge!(g,1,2)
addedge!(g,1,4; directed=false)
addedge!(g,2,3)
addedge!(g,4,2)
addedge!(g,4,5; directed=false)
addedge!(g,4,7; directed=false)
addedge!(g,5,6)
addedge!(g,7,6)
addedge!(g,7,5)
addedge!(g,8,3)

PDAGtoDAG!(g)
graphVStructure!(g)

meekRules!(g)

@test isneighbor(g,2,3) == true
@test isneighbor(g,8,3) == true


collect(edges(g))