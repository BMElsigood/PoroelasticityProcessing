module PoroelasticityProcessing

using DataFrames, Roots, SpecialFunctions, QuadGK, NumericalIntegration
using Statistics: mean
include("types.jl")
include("functions.jl")
include("diffusionrampsource.jl")
export keymechparams
export linBz, linB, linEz, volStrain, linK, linνz, linEν_u, darcyperm, gassmann, gassmannCijkl, axstresssteps,radstresssteps
export radstressstepsdiff,axstressstepsdiff,pressurerampsolution
end
