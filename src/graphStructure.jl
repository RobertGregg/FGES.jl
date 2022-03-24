####################################################################
# Basic Data Structures
####################################################################


#------PDAGs------#


#PDAG is parametric type because A could be a BitArray or a SubArray
"""
Define a Partially Directed Acyclic Graph (PDAG).

PDAGs are defined by the number of vertices and a Boolean adjacency matrix
"""
mutable struct PDAG{S<:AbstractMatrix{Bool}}
    nv::Int #number of vertices
    A::S #adjacency matrix to represent graph
end


"""
    PDAG(x::Integer)
Create an empty graph with `x` vertices.
"""
PDAG(nv::Integer) = PDAG(nv, falses(nv,nv))


"""
    PDAG(numNodes, numEdges; seed=-1)
Generate a random PDAG with `numNodes` vertices and `numEdges` edges.
"""
function PDAG(numNodes::Integer, numEdges::Integer; seed::Int=-1)

    #Calculate the maximum number of edges
    maxEdges = numNodes*(numNodes - 1) ÷ 2
    @assert(numEdges <= maxEdges, "Maximum number of edges for this graph is $maxe")

    #If no random seed is given, use the default RNG
    rng = seed==-1 ? Xoshiro() : Xoshiro(seed);

    #Generate a empty graph
    g = PDAG(numNodes)

    #Continually add edges until the correct number is reached
    edgesAdded = 0
    while edgesAdded < numEdges
        src = rand(rng, 1:numNodes)
        dest = rand(rng, 1:numNodes)
        orientation = rand(rng, Bool)
        
        #Do not allow self loops and do not over-write edges
        if src != dest && !isadjacent(g, src, dest)
            orientation ? addedge!(g, src, dest, directed=false) : addedge!(g, src, dest)
            edgesAdded += 1
        end
    end

    #Return the random PDAG
    return g
end


#------Edges------#

#A data type to hold edge Information
# parent → child 
"""
    Edge(parent, child, directed::Bool)
    Edge(Tuple(parent,child), directed::Bool)
    Edge(CartesianIndex(parent,child), directed::Bool)
A data type to hold edge information.

It is also possible to convert a `Tuple` or `CartesianIndex` into an Edge.
"""
struct Edge
    parent::Int
    child::Int
    directed::Bool
end

#converting edges from other types
Edge(idx::Tuple{Int, Int}, directed) = Edge(idx...,directed)
Edge(idx::CartesianIndex{2}, directed) = Edge(Tuple(idx),directed)

"A data type used to iterate through the edges of the graph"
struct IterEdge{G<:PDAG}
    g::G 
end



####################################################################
# Overload Base functions
####################################################################

#Used to subset the graph
#nodes contains a list of nodes you want to keep in the graph
#no need to make a copy of the graph using views!
view(g::PDAG, nodes) = PDAG(length(nodes), view(g.A, nodes, nodes))


#An element of an IterEdge is an Edge (gives functions like collect() information about the element type)
eltype(::Type{IterEdge{PDAG{T}}}) where T = Edge


#knowing the number of edges will allow us to pre-allocate a correctly sized collection of edges
#ne() is a function defined below to get the number of edges
length(iterEdge::IterEdge) = ne(iterEdge.g)


#overload to print PDAG (requires nv() and ne() which are functions defined below)
function show(io::IO, g::PDAG)
    print(io, "{$(nv(g)), $(ne(g))} PDAG")
end


#overload to print Edge type
function show(io::IO, edge::Edge)
    if edge.directed
        print(io,"$(edge.parent) → $(edge.child)")
    else
        print(io,"$(edge.parent) - $(edge.child)")
    end
end

####################################################################
# Number of Vertices and Edges
####################################################################

"Return the number of vertices"
nv(g::PDAG) = g.nv

"Return the number of edges"
function ne(g::PDAG)
    A = g.A
    totalEdges = 0
    #basically count each edge twice because some edges appear only on upper/lower diagonal but undirected edges appear on both
    for i in eachindex(A)
        if A[i] && A'[i]  #undirected
            totalEdges += 1
        elseif A[i] #directed
            totalEdges += 2
        end
    end

    return totalEdges ÷ 2 #integer division
end

####################################################################
# is* functions
####################################################################

"""
    isneighbor(g::PDAG, x, y)
    isneighbor(g::PDAG, edge::Edge)
Test if `x` and `y` are connected by a undirected edge in the graph `g`.

Can also perform the same test given an `edge`.
"""
isneighbor(g::PDAG, x, y) = g.A[x,y] && g.A[y,x]
isneighbor(g::PDAG, edge::Edge) = isneighbor(g, edge.parent, edge.child)


"""
    isadjacent(g::PDAG, x, y)
    isadjacent(g::PDAG, edge::Edge)
Test if `x` and `y` are connected by a any edge in the graph `g`.

Can also perform the same test given an `edge`.
"""
isadjacent(g::PDAG, x, y) = g.A[x,y] || g.A[y,x]
isadjacent(g::PDAG, edge::Edge) = isadjacent(g, edge.parent, edge.child)


"""
    isparent(g::PDAG, x, y)
    isparent(g::PDAG, edge::Edge)
Test if `x` is a parent of `y` in the graph `g`.

Can also perform the same test given an `edge`.
"""
isparent(g::PDAG, x, y) = g.A[x,y] && !g.A[y,x]
isparent(g::PDAG, edge::Edge) = isparent(g, edge.parent, edge.child)


"""
    ischild(g::PDAG, x, y)
    ischild(g::PDAG, edge::Edge)
Test if `x` is a child of `y` in the graph `g`.

Can also perform the same test given an `edge`.
"""
ischild(g::PDAG, x, y) = !g.A[x,y] && g.A[y,x]
ischild(g::PDAG, edge::Edge) = ischild(g, edge.parent, edge.child)

"""
    issink(g::PDAG, x)
Test if `x` a sink node in the graph `g` (i.e. there are no arrows pointing out of `x`).
"""
function issink(g::PDAG, x)
    for y ∈ vertices(g)
        if isparent(g,x,y)
            return false
        end
    end
    return true
end

"""
    issource(g::PDAG, x)
Test if `x` a source node in the graph `g` (i.e. there are no arrows pointing into `x`).
"""
function issource(g::PDAG, x)
    for y ∈ vertices(g)
        if ischild(g,x,y)
            return false
        end
    end
    return true
end


"""
    isclique(g::PDAG, nodes)
Return `true` if all vertices in `nodes` are connected to each other in the graph `g`.
"""
function isclique(g::PDAG, nodes::Vector{T}) where T<:Integer

    for (v₁, v₂) in combinations(nodes,2)
        if !isadjacent(g, v₁, v₂)
            return false
        end
    end

    return true
end


#Do the nodesRemoved block all semi-directed paths between src and dest?
"""
    isblocked(g::PDAG, src, dest; nodesRemoved=nothing)
Return `true` if there is no semi-directed path between `src` and `dest` in the graph `g`.

Optionally, a set of vertices (`nodesRemoved`) can be removed from the graph before searching for a semi-directed path.

A semi-directed path between `src` and `dest` is a list of edges in `g` where every edge is either undirected or points toward `dest`. 

    src → x₁ - x₂ → dest ✓
    src → x₁ ← x₂ - dest ✖
"""
function isblocked(g::PDAG, src, dest; nodesRemoved=nothing)

    #use breath first search (in hopes that path is short)

    #Create a graph subset where we remove the nodes in nodesRemoved
    if !isnothing(nodesRemoved)
        nodeSubset = filter(x-> x∉nodesRemoved, vertices(g))
        g = view(g, nodeSubset)
        src = searchsortedfirst(nodeSubset,src)
        dest = searchsortedfirst(nodeSubset,dest)
    end

    #first check if scr or dest have no neighbors
    if countNeighbors(g, src)==0 || countNeighbors(g, dest)==0
        return true
    end

    #Keep track of all the nodes visited and a queue of nodes to visit
    visited = falses(nv(g))
    queue = [src]

    #loop over all nodes in queue until it is empty
    while !isempty(queue)
        #Get the next node from the queue
        currentNode = popfirst!(queue)
        
        #if this edge was already visited from an undirected path, skip it
        #needed because of undirected edges
        visited[currentNode] && continue
        
        #Get the undirected and outgoing nodes from the current node
        nextNodes = children(g,currentNode) ∪ neighbors_undirect(g,currentNode)
        
        #Check if destination node was found
        if dest ∈ nextNodes
            return false
        end
        
        #If not, append the undirected and outgoing nodes to the queue
        append!(queue, nextNodes)
        
        #Mark the node as visited 
        visited[currentNode] = true
    end

    return true
end

"Return true if `edge` is directed"
isdirected(edge::Edge) = edge.directed

####################################################################
# Add and Remove 
####################################################################

"""
    addedge!(g::PDAG, x, y; directed=true)
    addedge!(g::PDAG, edge::Edge)
Add the edge `x`→`y` to the graph `g`.

Can also perform the same addition given an `edge`.

By default, this function will direct the new edge from `x` to `y` unless the `directed` keyword is changed to `false`.
"""
function addedge!(g::PDAG, x, y; directed=true)
    g.A[x,y] = true
    g.A[y,x] = !directed
    return nothing
end

#If an edge type if given, call the previous function
addedge!(g::PDAG, edge::Edge) = addedge!(g::PDAG, edge.parent, edge.child, directed=edge.directed)

"""
    remedge!(g::PDAG, x, y)
    remedge!(g::PDAG, edge::Edge)
Remove any edge between `x` and `y` in the graph `g`.

Can also perform the same deletion given an `edge`.

No warning is given if `edge` is not in the graph.
"""
function remedge!(g::PDAG, x, y)
    g.A[x,y] = false
    g.A[y,x] = false
    return nothing
end

#remove edge (Edge type)
remedge!(g::PDAG, edge::Edge) = remedge!(g, edge.parent, edge.child)



####################################################################
# Neighborhood functions 
####################################################################

#Loop through all the vertices and perform some neighborhood test on the given vertex
#Output all vertices that pass the given test
function neighborsGeneral(isFunction::Function, g::PDAG, x)
    neighborList = Int64[]
    for vᵢ in vertices(g)
        #check depends on desired list of neighbors
        #e.g. isFunction == isParent when we're searching for vᵢ→x
        if isFunction(g,vᵢ,x)
            push!(neighborList,vᵢ)
        end
    end

    return neighborList
end

#All vertices connected to x
neighbors(g::PDAG, x) = neighborsGeneral(isadjacent, g, x)

#All vertices pointing to x
neighbors_in(g::PDAG, x) = neighborsGeneral(isparent, g, x)

#All vertices pointing away from x
neighbors_out(g::PDAG, x) = neighborsGeneral(ischild, g, x)

#All vertices connected to x by an undirected edge
neighbors_undirect(g::PDAG, x) = neighborsGeneral(isneighbor, g, x)

#these are just aliases for the functions above
parents(g::PDAG, x) = neighbors_in(g, x)
children(g::PDAG, x) = neighbors_out(g, x)


####################################################################
# Counting Functions
####################################################################

#Very similar to the neighborsGeneral function except we only keep a count of the valid nodes
function countGeneral(isFunction::Function, g::PDAG,x)
    counter=0
    for vᵢ in vertices(g)
        if isFunction(g,vᵢ,x)
            counter += 1
        end
    end
    return counter
end

#All vertices connected to x
countNeighbors(g::PDAG, x) = countGeneral(isadjacent, g, x)

#All vertices pointing to x
countNeighbors_in(g::PDAG, x) = countGeneral(isparent, g, x)

#All vertices pointing away from x
countNeighbors_out(g::PDAG, x) = countGeneral(ischild, g, x)

#All vertices connected to x by an undirected edge
countNeighbors_undirect(g::PDAG, x) = countGeneral(isneighbor, g, x)

#these are just aliases for the functions above
countParents(g::PDAG, x) = countNeighbors_in(g, x)
countChildren(g::PDAG, x) = countNeighbors_out(g, x)

####################################################################
# all* functions
####################################################################

#Get every undirected edge in the entire graph g
allundirected(g::PDAG) = Edge.(findall(tril(g.A .& g.A')),false)

#Get every directed edge in the entire graph g
alldirected(g::PDAG) = [Edge(i,true) for i in CartesianIndices(g.A) if g.A[i] && !g.A'[i]]

#combine and sort all directed and undirected edges
alledges(g::PDAG) = sort([allundirected(g); alldirected(g)], by=x->(x.parent, x.child))



####################################################################
# Iterators
####################################################################

#Given an initial (i,j) index, output the next non-zero index in adjacency matrix
function iterate(iterEdge::IterEdge, state=(1,1))
    g = iterEdge.g
    n = nv(g)
    (row,col) = state

    #loop starting at the current state (row,col)
    for i = row:n
        start = row==i ? (col+1) : (i+1)
        for j = start:n

            #set the next potential state
            state = (i,j)

            #Check for an edge
            if g.A[i,j] && g.A'[i,j] #undirected
                return (Edge(state,false), state)
            elseif g.A[i,j] #forward edge
                return (Edge(state,true), state)
            elseif g.A'[i,j] #backward edge
                return (Edge(reverse(state),true), state)
            end
        end
    end

    return nothing
end

"""
    edges(g::PDAG)

Return an iterator to generate all edges within the graph `g`. Use 
"""
edges(g::PDAG) = IterEdge(g)


"""
    vertices(g::PDAG)

Return an iterator to generate all vertices within the graph `g`. 
"""
vertices(g::PDAG) = Base.OneTo(nv(g))


####################################################################
# Misc
####################################################################

"""
    orientedge!(g::PDAG, x, y)

Update the edge `x`-`y` to `x`→`y` in the graph `g`. 
"""
function orientedge!(g::PDAG, x, y)
    g.A[y,x] = false
    return nothing
end

# for x-y, get undirected neighbors of y and any neighbor of x
calcNAyx(g::PDAG, y::Integer, x::Integer) = setdiff(neighbors_undirect(g,y) ∩ neighbors(g,x), x)
calcNAyx(g::PDAG, edge::Edge) = calcNAyx(g, edge.child, edge.parent)

#for x-y, undirected neighbors of y not connected to x
calcT(g::PDAG, y::Integer, x::Integer) = setdiff(neighbors_undirect(g,y), neighbors(g,x), x)
calcT(g::PDAG, edge::Edge) = calcT(g, edge.child, edge.parent)

