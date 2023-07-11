#Using a directed graph that tracks forward and reverse edges from Graphs.jl

####################################################################
# Iterators
####################################################################

"""
    allpairs(v)
Iterate through all pairs of items in `v`
"""
allpairs(v) = Iterators.filter(i -> isless(i...), Iterators.product(v,v))


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
    isoriented(g, edge::Edge)
    isoriented(g, x, y)
Test if `x` and `y` are connected by a directed edge in the graph `g`, either x←y OR x→y.
Can also perform the same test given an `edge`.
"""
isoriented(g, edge) = has_edge(g,edge) ⊻ has_edge(g,reverse(edge)) # xor
isoriented(g, x, y) = has_edge(g,x,y) ⊻ has_edge(g,y,x)



"""
    isclique(g, nodes)
Return `true` if all vertices in `nodes` are connected to each other in the graph `g`.
"""
function isclique(g, nodes)

    for (x,y) in allpairs(nodes)
        if !isadjacent(g, x, y)
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
    #src needs to have a child
    isempty(children(g, src)) && return true
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
    add_edge!(g,x,y)
    rem_edge!(g,y,x) #TODO Might not be needed
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