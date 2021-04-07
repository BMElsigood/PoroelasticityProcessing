module PoroelasticityProcessing

using DataFrames
using Statistics: mean
include("types.jl")
include("functions.jl")
export keymechparams
export linBz, linB, linEz, volStrain, linK, linνz, linEν_u, darcyperm, gassmann, gassmannCijkl, axstresssteps,radstresssteps

end
