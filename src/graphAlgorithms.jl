####################################################################
# Using Meek's Rules to Update PDAG
####################################################################

#Revert a graph to undirected edges and unshielded colliders (i.e. parents not adjacent)
function graphVStructure!(g)
    
    #undirect an edge if it does not participate in an unshielded collider
    for edge in alledges(g)
        if edge.directed && allshielded(g,edge) 
            #undirect by adding reverse edge
            add_edge!(g, edge.child, edge.parent)
        end
    end
end

allshielded(g,x,y) = all(isadjacent(g, p, x) for p in parents(g, y) if p ≠ x)
allshielded(g,edge) = allshielded(g, edge.parent, edge.child)


function meekRules!(g)
    
    rulesFound = true

    while rulesFound
        
        rulesFound = false

        for edge in allUndirectedEdges(g)

            #For clarity extract the edge vertices
            (x, y) = edge.parent, edge.child

            #TODO Can multiple rules pass for the same edge?
            if R1(g,x,y) || R2(g,x,y) || R3(g,x,y)
                #Change x-y to x→y
                orientedge!(g, x, y)
                rulesFound = true
            elseif R1(g,y,x) || R2(g,y,x) || R3(g,y,x)
                #Change y-x to y→x
                orientedge!(g, y, x)
                rulesFound = true
            end
        end

    end

    return nothing
end


function R1(g, x, y)
    #given x-y, look for patterns that match v₁→x and not(v₁→y)
    for v₁ in parents(g,x)
        if !isadjacent(g,v₁,y)
            return true
        end
    end
    return false
end


function R2(g,x,y)
    #given x-y, look for patterns that match x→v₁→y
    for v₁ in children(g,x)
        if isparent(g,v₁,y)
            return true
        end
    end
    return false
end

function R3(g,x,y)
        
    #given x-y, find x-v₁→y and x-v₂→y and v₁-v₂
    for (v₁,v₂) in allpairs(neighbors(g,x))
        if isparent(g,v₁,y) && isparent(g,v₂,y) && !isadjacent(g,v₁,v₂)
            return true
        end
    end
    
    return false
end