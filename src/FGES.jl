module FGES
    
#Packages + Reason required
using LinearAlgebra #for subarrays (?)
using Statistics #mean
using Graphs #FGES works on DAGs
using Combinatorics #loop through combinations and permutations
using LRUCache #Caching scores
using Dates #record the current time
using DataFrames #To allow dataframe inputs to fges


#Files to include
include("graphStructure.jl")
include("graphAlgorithms.jl")
include("mainAlgorithm.jl")
include("graphIO.jl")

export 

#graphStructure.jl
    allpairs,
    isadjacent,
    isneighbor,
    isparent,
    ischild,
    isdescendent,
    isoriented,
    isclique,

    neighbors,
    parents,
    children,
    descendents,
    ancestors,

    isblocked,
    orientedge!,
    calculateNAyx,
    calculateTyx,

#graphAlgorithms.jl
    graphVStructure!,
    categorizeNeighbors,
    setCategory!,
    meekRules!,

#mainAlgorithm.jl
    CurrentState,
    fges,
    Insert!,
    Delete!,
    forwardSearch!,
    backwardSearch!,
    nextInsertEquivClass!,
    nextDeleteEquivClass!,
    findBestInsert!,
    findBestDelete!,
    score,
    calculateMSE,

#graphIO
    saveGraph,
    loadGraph,
    edgetable
end