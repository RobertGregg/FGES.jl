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
