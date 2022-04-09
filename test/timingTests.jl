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
@btime issink(g,x)
@btime issource(g,x)
@btime isclique($g,$vecNodes)
@btime isdirected(edge)
@btime no_path($g,$x,$y,$vecNodes)

@btime addedge!(g,edge)
@btime remedge!(g,edge)

@btime neighbors(g,x);
@btime neighbors_in(g,x);
@btime neighbors_out(g,x);
@btime neighbors_undirect(g,x);
@btime parents(g,x);
@btime children(g,x);

@btime countNeighbors(g,x);
@btime countNeighbors_in(g,x);
@btime countNeighbors_out(g,x);
@btime countNeighbors_undirect(g,x);
@btime countParents(g,x);
@btime countChildren(g,x);

@btime edges(g)
@btime vertices(g)

@btime orientedge!(g,x,y)
@btime calcNAyx($g,$x,$y)
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


@benchmark calcNAyx2($g,$x,$y)


function f(g,x,y,nodes)
    for i=1:1000
        no_path(g,x,y,nodes)
    end

    return nothing
end


@profview fges(data)

@time fges(data)



a = rand(1:1000,100)
b = rand(1:1000,100)

@btime Set(a)
b1 = Set(b)

@benchmark union($a,$b)
@benchmark union($a1,$b1)


using Combinatorics, BenchmarkTools

#Some complicated function that outputs a single Float64 score
f(x,y) = sin(x) + cos(y)

#Struct to hold the best score and arguments
Base.@kwdef mutable struct output
    x = 0
    y = 0
    score = 0.0
end

#function to loop over different valid arguments for f
function findBest(out,v)
    for (x,y) in combinations(v,2)
        newscore = f(x,y)
        if newscore > out.score
            out.x = x
            out.y = y
            out.score = newscore
        end
    end
end


#Example iteration
v = 1:50
out = output()
@benchmark findBest($out,$v)

h1(v) = maximum((f(x,y),x,y) for (x,y) in combinations(v,2))

@benchmark h($v)

#Parallel version 
function findBest_par(out,v)
    @sync begin
        for (x,y) in combinations(v,2)
            Threads.@spawn begin 
                newscore = f(x,y)
                if newscore > out.score
                    out.x = x
                    out.y = y
                    out.score = newscore
                end
            end
        end
    end
end



v = 1:10
out = output()
@btime findBest_par($out,$v)




using IterTools, BenchmarkTools, ThreadsX

f(x,y) = sin(x) + cos(y)

v = 1:100

h1(v) = maximum((f(x,y),x,y) for (x,y) in subsets(v,2))
h1(v)


h2(v) = ThreadsX.maximum((f(x,y),x,y) for (x,y) in Iterators.product(v,v))
h2(v)

h3(v) = ThreadsX.maximum((f(x,y),x,y) for (x,y) in Iterators.product(v,v) if x<y)
h3(v)

h4(v) = ThreadsX.maximum((f(x,y),x,y) for (x,y) in Iterators.filter(i -> isless(i...), Iterators.product(v,v)))
h4(v)

@benchmark h1($v)
@benchmark h2($v)
@benchmark h3($v)
@benchmark h4($v)


a =  (1,2)
isless(a...)


allPairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))



function myFindmax(v)
    @floop for (x,y) in allPairs(v)
        score = f(x,y)
        @reduce() do (bestx = -1; x), (besty = -1; y), (bestScore = 0.0; score)
            if bestScore < score
                bestx = x
                besty = y
                bestScore = score
            end
        end
    end

    return (bestScore,bestx,besty)
end



@benchmark h1($v)
@benchmark h3($v)
@benchmark myFindmax($v)

