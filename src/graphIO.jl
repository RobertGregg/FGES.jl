####################################################################
# Save/load graphs
####################################################################

"""
    saveGraph(fileName, g)

Save the output of fges to a text file. 
"""
function saveGraph(fileName, g, featureNames)
    open(fileName,"w") do io

        for (i,name) in enumerate(featureNames)
            println(io,i,":",name)
        end

        println(io,"")

        doubleCountedEdges = Vector{Tuple{Int64, Int64}}()

        for edge in edges(g)
            
            x, y = edge.src, edge.dst

            if (y,x) ∈ doubleCountedEdges
                continue
            end

            if isoriented(g, edge)
                println(io,"$x → $y")
            else
                println(io,"$x - $y")
                push!(doubleCountedEdges,(x, y))
            end

        end

    end
end

saveGraph(fileName, g) = saveGraph(fileName, g, vertices(g))


"""
    loadGraph(fileName)

Load a saved output of fges from a text file. 
"""
function loadGraph(fileName)

    #Let's check if the node names are numbers or strings
    parseOutput = tryparse(Int,last(split(first(eachline(fileName)),":")))

    if isnothing(parseOutput)
        nodeNames = String[]
    else
        nodeNames = Int[]
    end

    for line in eachline(fileName)
        if isempty(line)
            break
        end

        if isnothing(parseOutput)
            push!(nodeNames, String(last(split(line,":"))))
        else
            push!(nodeNames, parse(Int,last(split(line,":"))))
        end
    end


    #find the largest number in the file
    numFeatures = length(nodeNames)

    g = SimpleDiGraph(numFeatures)

    for line in eachline(fileName)

        sublines = split(line," ")
        #Check if line was split
        if length(sublines) == 1
            continue
        end

        if isnothing(parseOutput)
            v₁ = findfirst(isequal(sublines[1]), nodeNames)
            v₂ = findfirst(isequal(sublines[3]), nodeNames)
        else
            v₁ = parse(Int,sublines[1])
            v₂ = parse(Int,sublines[3])
        end
        
        add_edge!(g, v₁, v₂)

        #Check if edge is directed
        if sublines[2] ≠ "→"
            add_edge!(g, v₂, v₁)
        end
    end

    return g, nodeNames
end

####################################################################
# Convert graph to table
####################################################################

#Well that was easy, thanks Tables.jl
edgetable(g::T) where T<:Graphs.AbstractGraph = DataFrame(edges(g))