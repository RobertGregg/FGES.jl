using Revise
using FGES
using Random

Random.seed!(314)
data = rand(100,100)

@time fges(data)