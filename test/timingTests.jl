using Revise
using FGES, BenchmarkTools
using Random


g = PDAG(1000,10000)
x = 1
y = 2
vecNodes = [3,4,5]
edge = Edge(x,y,true)

@btime isneighbor(g,x,y)
@btime isadjacent(g,x,y)
@btime isparent(g,x,y)
@btime ischild(g,x,y)
@btime isclique($g,$vecNodes)
@btime isdirected(edge)
@btime isblocked($g,$x,$y,$vecNodes) #61.167 μs (6 allocations: 8.61 KiB)

@btime addedge!(g,edge)
@btime remedge!(g,edge)

@btime neighbors($g,$x);
@btime neighbors_in($g,$x);
@btime neighbors_out($g,$x);
@btime neighbors_undirect($g,$x);
@btime parents($g,$x);
@btime children($g,$x);
@btime descendent($g,$x);

@btime countNeighbors($g,$x);
@btime countNeighbors_in($g,$x);
@btime countNeighbors_out($g,$x);
@btime countNeighbors_undirect($g,$x);
@btime countParents($g,$x);
@btime countChildren($g,$x);

@btime edges($g);
@btime vertices($g);

@btime orientedge!($g,$x,$y);
@btime calcNAyx($g,$x,$y);
@btime degreeAverage,


#Create some  synthetic data for testing
Random.seed!(314)
numFeatures = 100
numObservations = 1000
data = zeros(numObservations, numFeatures)

for i in 1:numFeatures
    if i ≤ 5
        data[:,i] = randn(numObservations)
    else
        data[:,i] = sum(rand()*data[:,i-j] for j∈1:5) + randn(numObservations)
        #Stop the data from exploding
        data[:,i] ./= mean(data[:,i])
    end
end

@btime ParseData($data, nothing, 1.0) #6.614 μs (21 allocations: 39.33 KiB)


#Score Function
dataParsed = ParseData(data, nothing, 1.0)
#No Parents (calls mean)
nodeParents = Int[]; node=1
@btime score($dataParsed, $nodeParents, $node, false)

#Simple regression
nodeParents = 2; node=1
@btime score($dataParsed, $nodeParents, $node, false)

#Linear Solve (qr)
nodeParents = [2,3]; node=[1]
@btime score($dataParsed, $nodeParents, $node, false)

#Vectors of 1 allocate? no, but nearly double time
nodeParents = [2]; node=[1]
@btime score($dataParsed, $nodeParents, $node, false)


#graphAlgorithms.jl
@btime graphVStructure!($g);
@btime meekRules!($g); #17.958 μs (128 allocations: 16.50 KiB)

#mainAlgorithm.jl
fges,
ParseData,
Search!,
Insert!,
Delete!,
statusUpdate



function calcNAyx2(g::PDAG, y, x)
    neighborList = Int64[]

    #loop through all vertices except for x
    for vᵢ in Iterators.filter(!isequal(x),vertices(g))
        #look for neighbors of y and adjacencies to x
        if isneighbor(g,vᵢ,y) && isadjacent(g,vᵢ,x)
            push!(neighborList,vᵢ)
        end
    end

    return neighborList
end
