#Library to connect julia and R
library(JuliaConnectoR)

##################################
#Example 1: raw data matrix
##################################

#Load the FGES package in from julia
FGES <- juliaImport("FGES")
Base <- juliaImport("Base") #Base julia functions

#Create a random data matrix in R
m <- matrix(runif(10000*100) , ncol = 100)

#Run the fges algorithm in R from Julia
g <- FGES$fges(m, verbose=TRUE)

#Display the edges found by fges
#Julia returns lazy iterator so we need to collect the output into a vector
Base$collect( FGES$alledges(g) )

#Save the graph as a text file
FGES$saveGraph("myGraph.txt", g)

#Load the graph from a text file
gload <- FGES$loadGraph("myGraph.txt")

#Convert the edge list to an R dataframe
df <- as.data.frame(FGES$edgetable(gload))

##################################
#Example 2: Dataframes
##################################

#Load the FGES package in from julia
FGES <- juliaImport("FGES")

#Import DataFrames from Julia
DataFrames <- juliaImport("DataFrames")

#Create a random dataframe in R
df <- as.data.frame(matrix(runif(10000*100) , ncol = 100))

#Convert the R dataframe to a Julia DataFrame
df <- DataFrames$DataFrame(df)

#Pass the dataframe to FGES
g <- FGES$fges(df)

#Display the edges found by fges
Base$collect( FGES$alledges(g) )
