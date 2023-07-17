using Revise
using Graphs
using FGES, BenchmarkTools
using Random

Random.seed!(314)
g = SimpleDiGraph(1000,10000)
x = 1
y = 2
vecNodes = [3,4,5]
edge = Edge(x,y)

@btime allpairs($vecNodes);

@btime isneighbor($g,$x,$y)
@btime isadjacent($g,$x,$y)
@btime isparent($g,$x,$y)
@btime ischild($g,$x,$y)
@btime isclique($g,$vecNodes)
@btime isoriented($g,$edge)
@btime isblocked($g,$x,$y,$vecNodes) #allocates


@btime FGES.neighbors($g,$x);
@btime parents($g,$x);
@btime children($g,$x);
@btime descendents($g,$x);
@btime ancestors($g,$x);


@btime orientedge!($g,$x,$y);
@btime calculateNAyx($g,$x,$y);
@btime calculateTyx($g,$x,$y);

fges(rand(10,10))

