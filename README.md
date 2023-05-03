# Fast Greedy Equivalence Search (FGES)

Given a set of observational data, this algorithm will construct a directed acyclic graph where nodes represent variables in the data and edges represent causal connections between those variables.

## Installing the package

The package is still undergoing major changes, however it can be tested using 

```julia
pkg> add https://github.com/RobertGregg/FGES.jl.git
```

## Running the package

Given a matrix of data where each column is a variable and each row is an observation:

```julia
using FGES

numFeatures = 10
numObservations = 1000

data = rand(numObservations, numFeatures)
g = fges(data)
```

The `fges()` function can be run in debug mode, where additional information about the algorithm's progress is printed.

```julia
data = rand(50,50)
g = fges(data, debug=true)
```

You can also pre-compute a [scatter matrix](https://en.wikipedia.org/wiki/Scatter_matrix) for large systems. This is only recommended if this matrix will not fit into memory and you need to use [memory maps](https://docs.julialang.org/en/v1/stdlib/Mmap/). Data standardization is required.

```julia
using Statistics

data = rand(50,50)
#Copy the data and standardize each column
normAugData = copy(data)
normAugData .-= mean(normAugData, dims=1) #subtract mean
normAugData ./= std(normAugData, dims=1) #divide by standard deviation

#Augment a column of ones on the end for linear regression
normAugData = [normAugData ones(numObservations)]
scatterMat = normAugData * normAugData'
g = fges(data, scatterMat=scatterMat)
```

Finally, you can use a DataFrame as input and FGES will match column names to node labels
```julia
using DataFrames, Random

Random.seed!(951)
df = DataFrame(Height=rand(100), Age=rand(100), Weight=rand(100), HeartRate=rand(100))
g = fges(df)

collect(edges(g))
    3-element Vector{Edge{String}}:
    Height → Age
    Height - Weight
    HeartRate → Age
```

## Saving the Graph

The graph can saved as a text file using the following function:

```julia
data = rand(50,50)
g1 = fges(data)

fileName = "myGraph.txt"
saveGraph(fileName,g1)

#load that graph back in
g2 = loadGraph(fileName)
```

## Interfacing with R

This package can be used in R by connecting through Julia. This can be useful if you depend on packages exclusive to R. Note that you are still required to download Julia and the FGES package. 

```R
#Library to connect julia and R
library(JuliaConnectoR)

#Load the FGES package in from julia
FGES <- juliaImport("FGES")

#Create a random data matrix in R
m <- matrix(runif(10000*100) , ncol = 100)

#Run the fges algorithm in R from Julia
g <- FGES$fges(m)

#Display the edges found by fges
FGES$alledges(g)

#Save the graph as a text file
FGES$saveGraph("myGraph.txt", g)

#Load the graph from a text file
gload <- FGES$loadGraph("myGraph.txt")

#Convert the edge list to an R dataframe
df <- as.data.frame(FGES$edgetable(g))
```
