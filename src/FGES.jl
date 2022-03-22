module FGES

#Packages + Reason required
using LinearAlgebra #for subarrays (?)
using Statistics #mean, std
using Random #set seeds in generation of random graph
using Combinatorics #loop through combinations and permutations
using Memoization, LRUCache #Caching scores
using Revise

#Base functions to overload
import Base: show, view, iterate, eltype, length, == #functions from base to extend

end
