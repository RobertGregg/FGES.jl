####################################################################
# Graph Algorithms
####################################################################


function topological_sort(g::PDAG)
    sortedNodes, visited = Vector{Int}(), falses(g.nv)

    for node in vertices(g)
        if !visited[node]
            dfs(g,node,visited,sortedNodes)
        end
    end
    
    return sortedNodes
end


function dfs(g, currentNode, visited, sortedNodes)

    visited[currentNode] = true

    for parent in parents(g, currentNode)
        if !visited[parent]
            dfs(g,parent,visited,sortedNodes)
        end
    end

    push!(sortedNodes,currentNode)
end


function PDAGtoDAG!(g::PDAG)

    #From: https://github.com/juangamella/ges
    # NOTE!!!: Algorithm (1) is from the 1992 paper "A simple
    # algorithm to construct a consistent extension of a partially
    # oriented graph" by Dorit Dor and Michael Tarsi. There is an
    # ERROR in the summarized version in Chickering's paper. In
    # particular, the condition that N_x U Pa_x is a clique is not
    # equivalent to the condition from Dor & Torsi that every neighbor
    # of X should be adjacent to all of X's adjacent nodes. The
    # condition summarized in Chickering is more restrictive (i.e. it
    # also asks that the parents of X are adjacent to each other), but
    # this only results in an error for some graphs, and was only
    # uncovered during exhaustive testing.


    vertRemaining = collect(vertices(g))
    gSub = view(g,vertRemaining)
    #Instead of making a copy of the graph, create a container to save all the new edge orientations, then apply them at the end
    #gOut = deepcopy(g)
    edgesToApply = Vector{Edge}()

    while nv(gSub) > 1
        #Look for a sink node
        for x in vertices(gSub)
            if issink(gSub,x) 
                #caluclate it's undirected neighbors
                undirectNeighbors = neighbors_undirect(gSub,x)
                #Check if the sink neighbors are valid
                if validSink(gSub,x,undirectNeighbors)
                    #Direct all undirected neigbors into discovered edge
                    for y in undirectNeighbors
                        push!(edgesToApply, Edge(vertRemaining[y],vertRemaining[x],true))
                    end
                    #Remove the sink node from the graph and repeat
                    deleteat!(vertRemaining,x)
                    gSub = view(g,vertRemaining)
                    break
                end
            end
        end
    end

    #Apply all the changes found (addedge! overwites edges)
    for edge in edgesToApply
        addedge!(g, edge)
    end

    return g
end

#Test for all y-x that y is adjacent to all other vertices adjacent to x
function validSink(gSub,x, undirectNeighbors)

    #all nodes connected by undirected edges to x should form a clique
    if !isclique(gSub, undirectNeighbors)
        return false
    end

    #Every node with an undirected edge is connected to every parent of x
    for p in parents(gSub,x)
        for undir in undirectNeighbors
            if !isadjacent(gSub, p, undir)
                return false
            end
        end
    end

    return true
end

#After sorting the vertices, sort the edges 
function orderEdges(g::PDAG)

    sortedNodes = topological_sort(g)
    sortedEdges = Vector{Edge}()
    
    for node in topological_sort(g)
        currentParents = parents(g,node)
        #If the next sorted node has parents
        if !isempty(currentParents)
            #find how the parents are ordered
            sortedParents = filter(x->x∈currentParents,sortedNodes)
            #append the parents edges in reverse order
            append!(sortedEdges,Edge.(reverse(sortedParents), node, true))
        end
    end

    return sortedEdges
end

#Loop through sorted edges and label them
function labelEdges(g::PDAG)
    sortedEdges = orderEdges(g)
    labels = [:unknown for _ in eachindex(sortedEdges)]

    #Stop when all edges are labelled
    while any(labels .== :unknown)
        #Get the lowest ordered edge (x→y) with an unknown label
        currentEdge = sortedEdges[findfirst(isequal(:unknown),labels)]
        #This sub-routine was created to avoid a goto statement
        labelEdges_helper!(g, sortedEdges, currentEdge, labels)
    end

    return (sortedEdges,labels)
end

#Using a subroutine avoids a @goto statement in the original algorithm
function labelEdges_helper!(g, sortedEdges, currentEdge, labels)
     
     #loop through the sorted compelled edges
     for (i,thisEdge) in enumerate(sortedEdges)
         # try to find a compelled edge (w→x) where (w→x→y)
         if (labels[i] == :compelled) && (thisEdge.child == currentEdge.parent)
             #if w is not a parent of y...
             if !isparent(g, thisEdge.parent, currentEdge.child)
                 #Label x→y and every edge incident into y with “compelled” and exit
                 currentParents = findall(x-> x.child == currentEdge.child, sortedEdges)
                 labels[currentParents] .= :compelled
                 return nothing
             else
                #label w→y as compelled and continue search
                 labels[i] = :compelled
             end
         end
     end

     #get all edges where y is the child (currentEdge = x→y)
     currentParents = findall(x-> x.child == currentEdge.child, sortedEdges)
     #Compel currentEdge and all currentParents
     if length(currentParents) > 1
         labels[currentParents] .= :compelled
     #If no incident edges found, the currentEdge is reversible
     else
         labels[currentParents] .= :reversible
     end

     return nothing
end

"""
    PDAGtoCompletePDAG!(g::PDAG)
Convert a PDAG into a completed PDAG.

A complete PDAG has the nice property that any DAG extended from it will  be contained in the same equivalence class. This setp is required after every `Insert` or `Delete` in FGES.

The procedure runs through the following steps:

    PDAG → DAG → Labeled DAG → Complete PDAG
"""
function PDAGtoCompletePDAG!(g::PDAG)

    #Step 1: Convert the PDAG into a consistant DAG extension
    gDAG = PDAGtoDAG!(g)

    #Step 2: Label all the edges as either "compelled" or "reversible"
    (sortedEdges,labels) = labelEdges(gDAG)

    #Direct the compelled edges and undirect the reversible 
    for (edge, label) in zip(sortedEdges, labels)
        addedge!(g, edge.parent, edge.child, directed = label==:compelled ? true : false)
    end

    return g
end 



function meekRules(g::PDAG)
    
    #Loop through all the edges in the graph
    for edge in edges(g)
        if !edge.directed
            #Get the neighbors of the parent and child vertex
            parentNeighbors, childNeighbors =  get_categorized_neighbors(g, edge)


            #If both x and y have parents (p₁ → x - y ← p₂) do we just not orient or is the graph not extendable? 
        end
    end
end


function get_categorized_neighbors(g::PDAG, edge::Edge)
    x, y = edge.parent, edge.child

    parentNeighbors = Dict{Int, Symbol}()
    childNeighbors = Dict{Int, Symbol}()
    for vᵢ in vertices(g)
        if isadjacent(g, x, vᵢ)
            parentNeighbors[vᵢ] = set_category(g, x, vᵢ)
        end

        if isadjacent(g, y, vᵢ)
            childNeighbors[vᵢ] = set_category(g, y, vᵢ)
        end
    end

    return (parentNeighbors, childNeighbors)
end

function set_category(g::PDAG, x::Integer, y::Integer)
    if isparent(g,x,y)
        return :parent
    elseif ischild(g,x,y)
        return :child
    else
        return :undirected
    end
end