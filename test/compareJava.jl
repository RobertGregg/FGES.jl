# using FGES
using Statistics, Random
using CSV, DataFrames
using Combinatorics
using Graphs
using FGES

Random.seed!(3)

#Generate a dataset
numFeatures = 100
numObservations = 10000
data = zeros(numObservations, numFeatures)

#number of previous features to use in caluclating the next feature
n=5

for i in 1:numFeatures
    if i ≤ n
        data[:,i] = randn(numObservations)
    else
        data[:,i] = sum(rand()*data[:,i-j] for j∈1:n) + randn(numObservations)
        #Stop the data from exploding
        data[:,i] ./= mean(data[:,i])
    end
end

#Run the julia version of FGES and time it
gJulia = @time fges(data; verbose=true)


#Save the dataset to run for java
df = DataFrame(data, :auto)
CSV.write("test/testDatasets/features$(numFeatures)_observations$(numObservations).txt", df ,delim=" ")


#Read and parse the java results
gJava = SimpleDiGraph(numFeatures)
javaFilepath = "test/testDatasets/features$(numFeatures)_observations$(numObservations)_out.txt"

javaResult = filter(line -> occursin(r"x\d+ \D{3} x\d+",line),readlines(open(javaFilepath)))
javaResult = split.(javaResult," ")

javaEdges = [parse.(Int, filter.(isdigit, line[[2,4]])) for line in javaResult]

for (i,edge) in enumerate(javaEdges)
    add_edge!(gJava,edge...)

    if javaResult[i][3] ≠ "-->"
        add_edge!(gJava,reverse(edge)...)
    end
end


#Generate the true graph
gTrue = SimpleDiGraph(numFeatures)

for i=n+1:numFeatures
    for j=1:n
        add_edge!(gTrue,i-j, i)
    end
end

#Calculate precision and recall
function contingencyTable(g,gTrue)

    continTable = zeros(Int,2,2)

    for (x,y) in allpairs(vertices(gTrue))
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


continTableJulia = contingencyTable(gJulia,gTrue)
modelPrecision(continTableJulia)
modelRecall(continTableJulia)

continTableJava = contingencyTable(gJava,gTrue)
modelPrecision(continTableJava)
modelRecall(continTableJava)


#Hamming distance
sum(abs,adjacency_matrix(gJulia) .- adjacency_matrix(gTrue))
sum(abs,adjacency_matrix(gJava) .- adjacency_matrix(gTrue))