####################################################################
# Buffer to optimize X\y
####################################################################

#Saves a copy of R from QR decomposition to speedup future solves
mutable struct LinearSolveBuffer{T}
    R::Matrix{T}
    Rsub::Matrix{T}
    b::Vector{T}
    ŷ::Vector{T}

    function LinearSolveBuffer(data::AbstractMatrix{T}) where T
        
        n,m = size(data)
        _,R = qr(data)

        return new{T}(R, copy(R), zeros(T,m), zeros(T,n))
    end

end

####################################################################
# Data structures for caching scores
####################################################################

#Dictionary that saves scores based on what column subset were tested
#Has a max size and removes the least used entry if full
mutable struct SimpleLRUCache{K,V} <: AbstractDict{K,V}
    capacity::Int
    cache::OrderedDict{K,V}
end

#Constructor
function SimpleLRUCache{K,V}(n) where {K,V}
    d = OrderedDict{K,V}()
    sizehint!(d,n)
    return SimpleLRUCache(n,d)
end

Base.length(lru::SimpleLRUCache) = length(lru.cache)
Base.haskey(lru::SimpleLRUCache, key) = haskey(lru.cache,key)


function Base.getindex(lru::SimpleLRUCache, key)
    value = pop!(lru.cache, key)
    lru.cache[key] = value
    return value
end

function Base.setindex!(lru::SimpleLRUCache, value, key)
    if haskey(lru, key)
        pop!(lru.cache, key)
    end
    lru.cache[key] = value
    if length(lru.cache) > lru.capacity
        popfirst!(lru.cache)
    end
    return nothing
end

function Base.get!(default::Base.Callable, lru::SimpleLRUCache{K,V}, key::K) where {K,V}

    if haskey(lru,key)
        return lru[key]
    else
        value = default() #careful, no check here for type
        lru[key] = value
        return value
    end

end

#The iterator needs to run in reverse to put the most recent key/value first
function Base.iterate(lru::SimpleLRUCache)
    lru.cache.ndel > 0 && OrderedCollections.rehash!(lru.cache)
    n = length(lru)
    n < 1 && return nothing
    return (Pair(lru.cache.keys[n], lru.cache.vals[n]), n-1)
end

function Base.iterate(lru::SimpleLRUCache, i)
    i < 1 && return nothing
    return (Pair(lru.cache.keys[i], lru.cache.vals[i]), i-1)
end

#Enhance Base.show to include capacity
Base.summary(io::IO, lru::SimpleLRUCache{K,V}) where {K,V} = 
    print(io, "SimpleLRUCache{$(K), $(V)} with $(length(lru))/$(lru.capacity) entries")


####################################################################
# Data structure to track FGES Steps
####################################################################

#TODO maybe simplify by creating smaller structs like LinearSolveBuffer
mutable struct CurrentState{T}

    #data and properties
    data::Matrix{T}
    featureNames::Vector{String}
    numFeatures::Int
    numObservations::Int
    #R::Matrix{T} #TODO move from LinearSolveBuffer?

    penalty::T
    buffer::LinearSolveBuffer{T}

    #Information about the current step
    x::Int
    y::Int
    bestScore::T
    bestSubset::Vector{Int}
    scoreBoard::SimpleLRUCache{Tuple{Vector{Int}, Int}, T}
    stage::String
    verbose::Bool

    function CurrentState(data::AbstractMatrix{T}, featureNames, penalty, verbose) where T

        numObservations, numFeatures = size(data)

        #Append column of ones to data for intercept in linear regressions performed in scoring function
        data = [data ones(T,numObservations)]

        buffer = LinearSolveBuffer(data)

        #Default values
        #TODO could maybe use @kwdef
        x = 1
        y = 2
        bestScore = zero(T)
        bestSubset = zeros(Int,100) #guessing nodes will have >100 edges
        scoreBoard = SimpleLRUCache{Tuple{Vector{Int}, Int}, T}(5_000_000) #e.g. ([1,4],6) → 23.11
        stage = "Forward Search"

        return new{T}(
            data,
            featureNames,
            numFeatures,
            numObservations,
            penalty,
            buffer,
            x,
            y,
            bestScore,
            bestSubset,
            scoreBoard,
            stage,
            verbose
            )
    end
end

function Base.show(io::IO, state::CurrentState)

    printstyled(io,"Current State\n", bold=true, color=:blue)
    print(io, "Stage: ")
    printstyled(io, "$(state.stage)\n", color= state.stage == "Forward Search" ? :green : :red)
    println(io, "Score: $(state.bestScore)")
    println(io, "Edge: $(state.x)→$(state.y)")
    println(io, "Cache: $(round(100length(state.scoreBoard) / state.scoreBoard.capacity, digits=3))%")
    println(io,"----------------------------------")
end 