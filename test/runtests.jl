using Revise
using FGES
using Test


numFeatures = 50
numObservations = 10000
data = zeros(numObservations, numFeatures)

for i in 1:numFeatures
    if i ≤ 2
        data[:,i] = randn(numObservations)
    else
        data[:,i] = sum(rand()*data[:,i-j] for j∈1:2) + randn(numObservations)
    end
end

fges(data, debug=true)
