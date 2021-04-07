module PoroelasticityProcessing

using DataFrames
include("types.jl")
include("functions.jl")
export linBz, linB, linEz, volStrain, linK, linνz, linEν_u, darcyperm, gassmann, gassmannCijkl, axstresssteps,radstresssteps

end
