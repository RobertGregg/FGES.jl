using Revise
using FGES
using Test
using Random, Statistics
using DataFrames
using Graphs



@testset "Precision/Recall" begin
    Random.seed!(314)
    #Generate a dataset
    numFeatures = 20
    numObservations = 1000000
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

    #Generate the true graph
    gTrue = SimpleDiGraph(numFeatures)

    for i=n+1:numFeatures
        for j=1:n
            add_edge!(gTrue,i-j, i)
        end
    end

    g = fges(data)

    function contingencyTable(g,gTrue)

        continTable = zeros(Int,2,2)
    
        for (x,y) in allpairs(vertices(gTrue))
            trueEdge = isadjacent(gTrue,x,y)
            estEdge = isadjacent(g,x,y)
    
            if trueEdge & estEdge #true positive
                continTable[1,1] += 1
            elseif !trueEdge & !estEdge #true negative
                continTable[2,2] += 1
            elseif !trueEdge & estEdge #false positive
                continTable[2,1] += 1
            else #false negative
                continTable[1,2] += 1
            end
        end
    
        return continTable
    end

    modelPrecision(continTable) = continTable[1,1] / sum(continTable[:,1])
    modelRecall(continTable) = continTable[1,1] / sum(continTable[1,:])

    gtable = contingencyTable(g,gTrue)

    @test modelPrecision(gtable) > 0.75
    @test modelRecall(gtable) > 0.75
end


@testset "Save/Load Graph" begin
    fileName1 = "./test/testDatasets/savedGraph1.txt"
    fileName2 = "./test/testDatasets/savedGraph2.txt"
    Random.seed!(314)
    data = rand(100,100)

    #Testing raw data
    g1 = fges(data)

    saveGraph(fileName1,g1)

    g1Loaded, featureNames = loadGraph(fileName1)

    @test g1 == g1Loaded

    #Testing DataFrames as well
    g2, state2 = fges(DataFrame(data,:auto); returnState=true)

    saveGraph(fileName2, g2, state2.featureNames)

    g2Loaded, featureNames = loadGraph(fileName2)

    @test g2 == g2Loaded
end


@testset "Penality Values" begin
    Random.seed!(314)
    data = rand(100,100)

    g = fges(data)
    gpenalty = fges(data,penalty=2)

    @test ne(g) > ne(gpenalty)
end


@testset "DataFrames Test" begin
    Random.seed!(314)

    data = rand(100,20)

    g1 = fges(data)
    g2 = fges(DataFrame(data,:auto))

    @test ne(g1) == ne(g2)
end


@testset "Getting skeleton + v-structures" begin

    g1 = SimpleDiGraph(4)

    # 1 → 2 ← 3 → 4 
    add_edge!(g1,1,2)
    add_edge!(g1,3,2)
    add_edge!(g1,3,4)

    graphVStructure!(g1)

    @test isdirected(g1,3,4) == false

    #Same graph but connect 1-3 to remove collider
    g2 = SimpleDiGraph(4)
    add_edge!(g2,1,2)
    add_edge!(g2,3,2)
    add_edge!(g2,3,4)
    add_edge!(g2,1,3)
    add_edge!(g2,3,1)

    graphVStructure!(g2)

    @test all(!edge.directed for edge in alledges(g2)) == true

    #shielded and un-shielded colliders
    #1 → 2 ← 3 ∪ 1 - 3 - 4 → 2 
    g3 = SimpleDiGraph(4)
    add_edge!(g3,1,2)
    add_edge!(g3,3,2)
    add_edge!(g3,1,3)
    add_edge!(g3,3,1)
    add_edge!(g3,3,4)
    add_edge!(g3,4,3)
    add_edge!(g3,4,2)

    graphVStructure!(g3)

    @test !isdirected(g1,2,3) == false
end
