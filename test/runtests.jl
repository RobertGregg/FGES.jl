using Revise
using FGES
using Test
using Random, Statistics




@testset "Edge Tests" begin
    Random.seed!(314)
    #Generate a dataset
    numFeatures = 100
    numObservations = 10000
    data = zeros(numObservations, numFeatures)

    #number of previous features to use in caluclating the next feature
    n=5

    for i in 1:numFeatures
        if i ≤ n
            data[:,i] = randn(numObservations)
        else
            data[:,i] = sum(rand()*data[:,i-j] for j∈1:n) + randn(numObservations)
            #Stop the data from exploding
            data[:,i] ./= mean(data[:,i])
        end
    end

    g = fges(data)

    @test ne(g) == 458
    @test mapreduce(isdirected, +, edges(g)) == 455
end


@testset "Save/Load Graph" begin
    fileName = "./test/testDatasets/savedGraph.txt"
    Random.seed!(314)
    data = rand(100,100)

    g = fges(data)

    saveGraph(fileName,g)

    gLoaded = loadGraph(fileName)

    @test g == gLoaded
end




@testset "Penality Values" begin
    fileName = "./test/testDatasets/savedGraph.txt"
    Random.seed!(314)
    data = rand(100,100)

    g = fges(data)
    gpenalty = fges(data,penalty=2)

    @test ne(g) == 230
    @test ne(gpenalty) == 23
end