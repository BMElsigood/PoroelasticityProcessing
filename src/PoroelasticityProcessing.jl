module PoroelasticityProcessing

using DataFrames, Roots, SpecialFunctions, QuadGK, NumericalIntegration, NumericalDifferentiation, VTIModuli
using Statistics: mean
include("types.jl")
include("functions.jl")
include("diffusionrampsource.jl")
include("diffusionplotting.jl")
export plotcheckfit,differrors
include("velocities.jl")
export keymechparams
export lininterp,linfit,linBz, linB, linEz, volStrain, linK, linνz, linEν_u, darcyperm,density, gassmann, gassmannCijkl, axstresssteps,radstresssteps,numdiff,alllinfit
export radstressstepsdiff,axstressstepsdiff,pressurerampsolution,psolvect
export svyCompliancet, surveytocompliancedf, Bx, Bz, B, sayers95, BxBzWong
end
