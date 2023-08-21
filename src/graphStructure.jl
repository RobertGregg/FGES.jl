#Using a directed graph that tracks forward and reverse edges from Graphs.jl

####################################################################
# Custom Edges
####################################################################

#The Graphs.jl edges() function outputs 1→2 and 2→1 if edge is undirected
#Here we're making a new iterator to avoid that
#In theory, I shoul define a new graph type, but that would be involved

"""
    GraphEdge(parent, child, directed::Bool)
    GraphEdge(Tuple(parent,child), directed::Bool)
    GraphEdge(CartesianIndex(parent,child), directed::Bool)
A data type to hold edge information.
"""
struct GraphEdge
    parent::Int
    child::Int
    directed::Bool
end


"A data type used to iterate through the edges of the graph"
struct IterEdge{G<:Graphs.AbstractGraph}
    g::G 
end

#(Remove the double counted edges)
Base.length(iterEdge::IterEdge) = ne(iterEdge.g) - mapreduce(x->!isdirected(iterEdge.g,x), +, edges(iterEdge.g)) ÷ 2

Base.eltype(::Type{IterEdge{T}}) where T = GraphEdge

#overload to print Edge type
function Base.show(io::IO, edge::GraphEdge)
    if edge.directed
        print(io,"$(edge.parent) → $(edge.child)")
    else
        print(io,"$(edge.parent) - $(edge.child)")
    end
end


####################################################################
# Iterators
####################################################################

"""
    allpairs(v)
Iterate through all pairs of items in `v`
"""
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))


#Given an initial (i,j) index, output the next non-zero index in adjacency matrix
function Base.iterate(iterEdge::IterEdge, state=(1,0))
    g = iterEdge.g
    (src, i) = state
    
    #Get the next edge skipping over unidrected edges
    while src ≤ nv(g)
        
        adjacencies = g.fadjlist[src]

        #If there are no adjacencies, go to next vertex
        if isempty(adjacencies)
            src += 1
            continue
        end

        #Increment to the next state
        if i ≥ length(adjacencies)
            i = 0
            src += 1
            continue
        else
            i += 1
        end

        dst = adjacencies[i]

        
        #We know src→dst, now check if dst→src
        undirectedEdge = insorted(src, g.fadjlist[dst])


        #Check if edge was already found
        if src > dst && undirectedEdge
            continue
        end
        
        return (GraphEdge(src, dst, !undirectedEdge), (src,i))
    end
    
    return nothing
end

"""
    alledges(g)

Return an iterator to generate all edges within the graph `g`. Similar to `Graphs.edges()` but does not double count undirected edges
"""
alledges(g) = IterEdge(g)
allUndirectedEdges(g) = Iterators.filter(edge->!edge.directed, alledges(g))

####################################################################
# Relationship between two verticies
####################################################################

"""
    isadjacent(g, x, y)
Test if `x` and `y` are connected by any edge in the graph `g`.
"""
isadjacent(g, x, y) = has_edge(g,x,y) || has_edge(g,y,x)


"""
    isneighbor(g, x, y)
Test if `x` and `y` are connected by a undirected edge in the graph `g`.
"""
isneighbor(g, x, y) = has_edge(g,x,y) && has_edge(g,y,x)


"""
    isparent(g, x, y)
Test if `x` is a parent of `y` in the graph `g`.
"""
isparent(g, x, y) = has_edge(g,x,y) && !has_edge(g,y,x)


"""
    ischild(g, x, y)
Test if `x` is a child of `y` in the graph `g`.
"""
ischild(g, x, y) = !has_edge(g,x,y) && has_edge(g,y,x)



"""
    isdescendent(g, x, y)
Return `true` if `x`←`y` OR `x`-`y` in the graph `g`.
"""
isdescendent(g, x, y) = has_edge(g,y,x)


"""
    isdirected(g, edge::Edge)
    isdirected(g, x, y)
Test if `x` and `y` are connected by a directed edge in the graph `g`, either x←y OR x→y.
Can also perform the same test given an `edge`.
"""
isdirected(g, edge) = has_edge(g,edge) ⊻ has_edge(g,reverse(edge)) # xor
isdirected(g, x, y) = has_edge(g,x,y) ⊻ has_edge(g,y,x)



"""
    isclique(g, nodes)
Return `true` if all vertices in `nodes` are undirected neighbors in the graph `g`.
"""
function isclique(g, nodes)

    for (x,y) in allpairs(nodes)
        if !isneighbor(g, x, y)
            return false
        end
    end

    return true
end


####################################################################
# Neighborhood functions 
####################################################################

#All undirected vertices to x
neighbors(g, x) = Iterators.filter(v -> v ∈ g.fadjlist[x], g.badjlist[x])

#All vertices pointing to x
parents(g, x) = Iterators.filter(v -> v ∉ g.fadjlist[x], g.badjlist[x])

#All vertices pointing away from x
children(g, x) = Iterators.filter(v -> v ∉ g.badjlist[x], g.fadjlist[x])

#All vertices that are connected but not a parent
descendents(g, x) = g.fadjlist[x]

#All vertices that are connected but not a child
ancestors(g, x) = g.badjlist[x]

#All vertices connected to x
adjacencies(g, x) = union(g.fadjlist[x], g.badjlist[x])


####################################################################
# FGES Specific functions
####################################################################

"""
    isblocked(g, src, dst, nodesRemoved)
Return `true` if there is no semi-directed path between `src` and `dst` in the graph `g`.

A set of vertices (`nodesRemoved`) can be removed from the graph before searching for a semi-directed path.

A semi-directed path between `src` and `dst` is a list of edges in `g` where every edge is either undirected or points toward `dst`. 

    src → x₁ - x₂ → dst ✓
    src → x₁ ← x₂ - dst ✖
"""
function isblocked(g, src, dst, nodesRemoved)

    #Keep track of all the nodes visited
    visited = falses(nv(g))

    # mark excluded vertices as visited
    for vᵢ in nodesRemoved 
        visited[vᵢ] = true
    end

    #For there to be a semi-directed path...
    #src needs to have a descendent not in nodesRemoved
    isempty(Iterators.filter(x->x∉nodesRemoved, descendents(g,src))) && return true
    #dst needs an ancestor not in nodesRemoved
    isempty(Iterators.filter(x->x∉nodesRemoved, ancestors(g,dst))) && return true

    queue = collect(children(g,src))
    visited[queue] .= true
    while !isempty(queue)
        currentNode = popfirst!(queue) # get new element from queue
        for vᵢ in descendents(g,currentNode)
            vᵢ == dst && return false
            if !visited[vᵢ]
                push!(queue, vᵢ) # push onto queue
                visited[vᵢ] = true
            end
        end
    end

    return true

end

"""
    orientedge!(g, x, y)
Update the edge `x`-`y` to `x`→`y` in the graph `g`. 
"""
function orientedge!(g, x, y)

    #Check if x-y exists
    !isneighbor(g,x,y) && return nothing

    #Remove the backwards edge to direct
    rem_edge!(g,y,x) 
    return nothing
end


"""
    calculateNAyx(g,y,x)
Verticies that are undirected neighbors of y and connect to x
"""
calculateNAyx(g,y,x) = Iterators.filter(vᵢ -> isadjacent(g,vᵢ,x), neighbors(g,y))


"""
    calculateTyx(g,y,x)
Verticies that are undirected neighbors of y and NOT connect to x
"""
calculateTyx(g,y,x) = Iterators.filter(vᵢ -> !isadjacent(g,vᵢ,x), neighbors(g,y))