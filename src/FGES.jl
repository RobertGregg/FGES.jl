module FGES
    
#Packages + Reason required
using LinearAlgebra #for qr
using Statistics #mean
using Graphs #FGES works on DAGs
using Combinatorics #loop through combinations and permutations
using OrderedCollections #Caching scores
using DataFrames #To allow dataframe inputs to fges
using ThreadsX #parallelization
using GLMakie, NetworkLayout #plotting PDAGs



#Files to include
include("dataStructures.jl")
include("graphStructure.jl")
include("graphAlgorithms.jl")
include("score.jl")
include("mainAlgorithm.jl")
include("graphIO.jl")
include("graphPlotter.jl")


#For now export everything to make testing easier
export 

#dataStructures.jl
    CurrentState,
    LinearSolveBuffer,

#graphStructure.jl
    allpairs,
    alledges,
    isadjacent,
    isneighbor,
    isparent,
    ischild,
    isdescendent,
    isdirected,
    isclique,

    neighbors,
    parents,
    children,
    descendents,
    ancestors,
    adjacencies,

    isblocked,
    orientedge!,
    calculateNAyx,
    calculateTyx,

#graphAlgorithms.jl
    graphVStructure!,
    categorizeNeighbors,
    setCategory!,
    meekRules!,

#score.jl
    score,
    calculateMSE,

#mainAlgorithm.jl
    fges,
    Insert!,
    Delete!,
    forwardSearch!,
    backwardSearch!,
    nextInsertEquivClass!,
    nextDeleteEquivClass!,
    findBestInsert!,
    findBestDelete!,

#graphIO
    saveGraph,
    loadGraph,
    edgetable,

#graphPlotter
    plotNetwork
end