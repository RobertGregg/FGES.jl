####################################################################
# Buffer to optimize X\y
####################################################################

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
# Data structures for algorithm
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
    scoreBoard::LRU{Tuple{Vector{Int}, Int}, T}
    stage::String
    verbose::Bool

    function CurrentState(data::AbstractMatrix{T}, featureNames, penalty, verbose) where T

        numObservations, numFeatures = size(data)

        #Append column of ones to data for intercept in linear regressions performed in scoring function
        data = [data ones(numObservations)]

        buffer = LinearSolveBuffer(data)

        #Default values
        #TODO could maybe use @kwdef
        x = 1
        y = 2
        bestScore = 0.0
        bestSubset = zeros(100) #guessing nodes will have >100 edges
        scoreBoard = LRU{Tuple{Vector{Int}, Int}, T}(maxsize=5_000_000) #e.g. ([1,4],6) → 23.11
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
    println(io, "Stage: $(state.stage)")
    println(io, "Score: $(state.bestScore)")
    println(io, "Edge: $(state.x)→$(state.y)")
    println(io, "Cache: $(round(100state.scoreBoard.currentsize / state.scoreBoard.maxsize, digits=3))%")
    println(io,"----------------------------------")
end 