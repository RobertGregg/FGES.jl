using Revise
using FGES, BenchmarkTools


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
@btime isblocked($g,$x,$y,$vecNodes) #allocates

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

#graphAlgorithms.jl
graphVStructure!,
categorizeNeighbors,
setCategory!,
meekRules!,

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
