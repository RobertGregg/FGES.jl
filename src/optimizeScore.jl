using BenchmarkTools
using Random
Random.seed!(31)

####################################################################
# Scoring function
####################################################################

#The goal here is to improve the scoring function by reducing allocations, parallelizition, and pruning

#Data to be used for testing
numObservations = 1000
numFeatures = 20

data = rand(numObservations,numFeatures)
data[:,end] .= 1.0

nodeParents = [1,5,8]
node = 10

function getRandNodes(numFeatures)
    n = rand(1:numFeatures-1)
    cols = randperm(numFeatures-1)

    nodeParents = sort(cols[1:n-1])
    node = cols[end]
    return nodeParents, node
end

#Here is a baseline scoring function to compare performance
function score(data, nodeParents, node)
    
    #Calculate datasize
    n, numFeatures = size(data)

    push!(nodeParents, numFeatures)

    @views begin
        X = data[:,nodeParents]
        y = data[:,node]
    end

    b = X\y

    ŷ = X*b

    
    mse = sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n

    k = length(nodeParents)
    p=1.0

    pop!(nodeParents)

    return -p*k*log(n) - n*log(mse)
end


####################################################################
# Benchmark and profiling
####################################################################

score(data,nodeParents,node)

@benchmark score($data,$nodeParents,$node)


function scoreLoop(data)
    for _ in 1:1000
        nodeParents,node = getRandNodes(20)
        score(data,nodeParents,node)
    end
    
    return nothing
end

@benchmark scoreLoop($data)

@profview scoreLoop(data)


####################################################################
# score with buffer for b and ŷ
####################################################################
using LinearAlgebra

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

#Loop through lower triangular matrix
#the 1st axis is reversed because we're moving non-zeros *up* above diagonal
lowTriIter(ax1, ax2) = Iterators.filter(i -> first(i)>last(i), Iterators.product(Iterators.reverse(ax1),ax2))
lowTriIter(A::AbstractMatrix) = lowTriIter(axes(A)...)


function qless!(b,X,y,R,Rsub,columnsKeep)

    numCol = length(columnsKeep)

    #Rsub = copy(R[:,columnsKeep]) #dont want to destroy R
    #copyto!(Rsub, R[:,columnsKeep])
    for (k,j) in enumerate(columnsKeep)
        @views Rsub[:,k] = R[:,j]
    end

    #Use a givens rotation to zero out values below diagonal
    @inbounds for (i,j) in lowTriIter(Rsub)        
        if Rsub[i,j] ≠ 0.0
            G,r = givens(Rsub,i-1,i,j)
            Rsub[i-1,j] = r
            Rsub[i,j] = 0.0

            for k in j+1:numCol
                (r1, r2) = Rsub[i-1,k], Rsub[i,k]
                Rsub[i-1,k] = G.c*r1 + G.s*r2
                Rsub[i,k] = -G.s*r1 + G.c*r2
            end
        end
    end

    Rv = UpperTriangular(view(Rsub,1:numCol,1:numCol))

    #Solve R'R*b = X'y
    mul!(b,X',y)
    ldiv!(Rv',b)
    ldiv!(Rv,b)

    return nothing
end

function score(data, nodeParents, node, buffer::LinearSolveBuffer)
    
    #Calculate datasize
    n, numFeatures = size(data)
    
    #Penalty for BIC
    p=1.0

    push!(nodeParents, numFeatures)
    k = length(nodeParents)

    @views begin
        X = data[:,nodeParents]
        y = data[:,node]
        b = buffer.b[1:length(nodeParents)]
        ŷ = buffer.ŷ
    end

    #b = X\y
    qless!(b, X, y, buffer.R, buffer.Rsub, nodeParents)
    
    #ŷ = X*b
    mul!(ŷ, X, b)

    
    mse = sum((yᵢ - ŷᵢ)^2 for (yᵢ, ŷᵢ) in zip(y,ŷ)) / n


    pop!(nodeParents)

    return -p*k*log(n) - n*log(mse)
end


####################################################################
# Benchmark and profiling
####################################################################

buffer = LinearSolveBuffer(data)
score(data,nodeParents,node,buffer)

@benchmark score($data,$nodeParents,$node,$buffer)

function scoreLoop(data,buffer)
    for _ in 1:1000
        nodeParents,node = getRandNodes(20)
        score(data,nodeParents,node,buffer)
    end

    return nothing
end

@benchmark scoreLoop($data,$buffer)

@profview scoreLoop(data,buffer)
