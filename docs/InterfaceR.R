#Library to connect julia and R
library(JuliaConnectoR)

#Load the FGES package in from julia
FGES <- juliaImport("FGES")

#Create a random data matrix in R
m <- matrix(runif(10000*100) , ncol = 100)

#Run the fges algorithm in R from Julia
g <- FGES$fges(m)

#Display the edges found by fges
FGES$alledges(g)

#Save the graph as a text file
FGES$saveGraph("myGraph.txt", g)

#Load the graph from a text file
gload <- FGES$loadGraph("myGraph.txt")

#Convert the edge list to an R dataframe
df <- as.data.frame(FGES$edgetable(g))
