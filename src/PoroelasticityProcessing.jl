module PoroelasticityProcessing

using DataFrames, Roots, SpecialFunctions, QuadGK, NumericalIntegration, VTIModuli
using Statistics: mean
include("types.jl")
include("functions.jl")
include("diffusionrampsource.jl")
include("velocities.jl")
export keymechparams
export linBz, linB, linEz, volStrain, linK, linνz, linEν_u, darcyperm,density, gassmann, gassmannCijkl, axstresssteps,radstresssteps
export radstressstepsdiff,axstressstepsdiff,pressurerampsolution
export svyCompliancet, surveytocompliancedf, Bx, Bz, B, sayers95
end
