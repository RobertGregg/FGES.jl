####################################################################
# Using Meek's Rules to Update PDAG
####################################################################

#Revert a graph to undirected edges and unshielded colliders (i.e. parents not adjacent)
function graphVStructure!(g::PDAG)
    
    #loop through all vertices
    for x in vertices(g)
        #get the parents of current vertex
        parentsX = parents(g,x)
        #if all the parents are adjacent, undirect the edges
        #(if there is only 1 parent, then it is still a clique)
        if isclique(g, parentsX)
            for p in parentsX
                addedge!(g, x, p, directed=false)
            end
        end
    end
end


#Thoughts on improvement
    #Would it be cleaner to find all unique triples where one edge is undirected?
    #Can this update be more local? Most edges will remain unchanged.

#This is set up so we only loop through all the vertices and edges once
function meekRules!(g::PDAG)
    
    #Loop through all the edges in the graph (x-y)
    for edge in edges(g)
        #We only need to update undirected edges
        if !edge.directed

            #For clarity extract the edge vertices
            x, y = edge.parent, edge.child

            #Label the neighbors of the nodes that comprise the current edge
            #Neighbors can be labelled "parent", "child", or "undirected"
            #e.g. v₁→x-v₂, v₁ is a parent and v₂ is undirected
            xNeighbors, yNeighbors =  categorizeNeighbors(g, x, y)

            #Rule 1: If x or y have a unique parent, direct edge away from it
            #Check if x has a unique parent 
            if R1(xNeighbors, yNeighbors)
                #Add x→y
                addedge!(g, x, y)
            #Check if y has a unique parent
            elseif R1(yNeighbors, xNeighbors)
                #Add y→x
                addedge!(g, y, x)
            end

            #Rule 2: Direct the edge away from a potential cycle
            #Check if x has a child that is a parent of y
            if R2(xNeighbors, yNeighbors)
                #Add x→y
                addedge!(g, x, y)
            #Check if y has a child that is a parent of x
            elseif R2(yNeighbors, xNeighbors)
                #Add y→x
                addedge!(g, y, x)
            end

            #Rule 3: Double Triangle, diagonal
            if R3(xNeighbors, yNeighbors)
                #Add x→y
                addedge!(g, x, y)
            elseif R3(yNeighbors, xNeighbors)
                #Add y→x
                addedge!(g, y, x)
            end

            #Rule 4: Double Triangle, side
            #Requires an additional test that can't be check with the neighbor dictionaries 
            #So, we need to pass the graph through as well
            if R4(xNeighbors, yNeighbors, g)
                #Add x→y
                addedge!(g, x, y)
            elseif R4(yNeighbors, xNeighbors, g)
                #Add y→x
                addedge!(g, y, x)
            end
            
        end
    end

    return nothing
end


function R1(neighborSet1, neighborSet2)
    #given x-y, look for patterns that match v₁→x and not(v₁→y)
    for (v₁, category) in neighborSet1
        if category == :parent && !haskey(neighborSet2,v₁)
            return true
        end
    end
    return false
end

function R2(neighborSet1, neighborSet2)
    #given x-y, look for patterns that match x→v₁→y
    for (v₁, category) in neighborSet1
        if category == :child && haskey(neighborSet2,v₁) && neighborSet2[v₁]==:parent
            return true
        end
    end
    return false
end

function R3(neighborSet1, neighborSet2)
        
    #given x-y, count the number of patterns that match x-v₁→y
    counter = 0
    for (v₁, category) in neighborSet1
        if category == :undirected && haskey(neighborSet2,v₁) && neighborSet2[v₁]==:parent
            counter += 1
            #If we find two such paths, the Rule 3 pattern is matched
            if counter == 2
                return true
            end
        end
    end

    return false
end

function R4(neighborSet1, neighborSet2, g)
        
    #given x-y, look for patterns that match x-v₁→v₂→y and x-v₂
    for (v₂, category) in neighborSet1
        #first look for x-v₂→y
        if category == :undirected && haskey(neighborSet2,v₂) && neighborSet2[v₂]==:parent
            #check for x-v₁ and v₁→v₂
            for (v₁, category) in neighborSet1
                if category == :undirected && isparent(g, v₁, v₂)
                    return true
                end
            end
        end 
    end

    return false
end


#Categorize neighboring edges as "parent", "child", or "undirected"
function categorizeNeighbors(g::PDAG, x, y)
    
    #Create a dictionary of neighbors (e.g. vertex 5 => :parent)
    xNeighbors = Dict{Int, Symbol}()
    yNeighbors = Dict{Int, Symbol}()

    #Loop through all the vertices in the graph
    for vᵢ in vertices(g)

        #If vᵢ is adjacent, give it a category
        if isadjacent(g, x, vᵢ)
            xNeighbors[vᵢ] = setCategory(g, vᵢ, x)
        end

        #Same procedure to the other vertex in the current edge
        if isadjacent(g, y, vᵢ)
            yNeighbors[vᵢ] = setCategory(g, vᵢ, y)
        end
    end

    return (xNeighbors, yNeighbors)
end


function setCategory(g::PDAG, v₁, v₂)

    #Test if v₁ → v₂
    if isparent(g,v₁,v₂)
        return :parent
    #Test if v₁ ← v₂
    elseif ischild(g,v₁,v₂)
        return :child
    #If not then we must have an undirected edge
    else
        return :undirected
    end
end

