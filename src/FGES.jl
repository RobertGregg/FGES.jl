module FGES
    
#Packages + Reason required
using LinearAlgebra #for subarrays (?)
using Statistics #mean, std
using Random #set seeds in generation of random graph
using Combinatorics #loop through combinations and permutations
using Memoization, LRUCache #Caching scores
#using Revise

#Base functions to overload
import Base: show, view, iterate, eltype, length, == #functions from base to extend

#Files to include
include("graphStructure.jl")
include("graphAlgorithms.jl")
include("mainAlgorithm.jl")

export 

#graphStructure.jl
    PDAG,
    Edge,

    nv,
    ne,

    isneighbor,
    isadjacent,
    isparent,
    ischild,
    issink,
    issource,
    isclique,
    isblocked,
    isdirected,

    addedge!,
    remedge!,

    neighbors,
    neighbors_in,
    neighbors_out,
    neighbors_undirect,
    parents,
    children,

    allundirected,
    alldirected,
    alledges,

    edges,
    vertices,

#graphStructure.jl
    PDAGtoCompletePDAG!,

#mainAlgorithm.jl
    fges


end