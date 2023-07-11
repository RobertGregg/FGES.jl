using GLMakie, NetworkLayout

#Smaller r will result in shorter edges
function nudge(r,p1,p2)
    θ =  atan(last(p2-p1), first(p2-p1))
    return p1 + [r*cos(θ), r*sin(θ)]
end


function plotNetwork(g; arrowgap=0.1, arrowsize=20, nodesize=30, layoutmethod=shell)
    
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1])
    hidedecorations!(ax)
    hidespines!(ax)

    
    #Determine the arrowheads
    heads =  fill('▲', ne(g))
    tails = [isoriented(g, edge) ? ' ' : '▲' for edge in edges(g)]
    
    #Calculate the node positions
    nodePositions = layoutmethod(adjacency_matrix(g))

    edgePairs = [(edge.src, edge.dst) for edge in edges(g)]
    start = [nudge(arrowgap, nodePositions[i], nodePositions[j]) for (i,j) in edgePairs]
    stop = [nudge(arrowgap, nodePositions[j], nodePositions[i]) for (i,j) in edgePairs]
    step = Vec2f.(stop .- start)

    #Plot the edges
    arrows!(ax, start, step; linecolor=(:black, 0.0), arrowhead=heads, arrowsize=arrowsize)
    arrows!(ax, stop, -step; arrowhead=tails, arrowsize=arrowsize)
    

    #Plot the nodes and add labels
    scatter!(ax, nodePositions, color=:dodgerblue, markersize=nodesize, strokewidth=1)
    text!(ax, nodePositions, text = string.(vertices(g)), align = (:center,:center))

    #Add some padding to the plot
    xlims!(ax,extrema(first.(nodePositions)) .* 1.2)
    ylims!(ax,extrema(last.(nodePositions)) .* 1.2)

    return fig
end