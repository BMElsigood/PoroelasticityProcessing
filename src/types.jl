#=e.g.

mutable struct AEparams
    folder::String
    masterfile::String
    sensorfile::String
    skipcomponent::Array{String,1}
    nchan::Int
    T0::Float64
    f0::Float64
    filter#::T where T<:FilterCoefficients
    fcc::Float64
    I_master::Int
    xcfront::Int
    xcback::Int
    Vmin::Float64
    Vmax::Float64
    pth::Paths
end

function show(io::IO, p::AEparams)
    println("AEparams:")
    for f in fieldnames(typeof(p))
        println(" .",f, " : ", getproperty(p,f))
    end
end

=#

"""
    choose mechdata we want for calc, could average sensors for pp and strain
"""
struct keymechparams
    time::Array{Float64,1}
    stress::Array{Float64,1}
    Pc::Array{Float64,1}
    pp::Array{Float64,1}
    εz::Array{Float64,1}
    εx::Array{Float64,1}
end
#=
function show(io::IO, data::keymechparams)
    println("keymechparams:")
    for f in fieldnames(typeof(data))
        println(" .",f, " : ", getproperty(data,f))
    end
end
=#
mutable struct undrainedaxstresssteps
    meanstress::Array{Float64,1}
    Bz::Array{Float64,1}
    Ez::Array{Float64,1}
    νz::Array{Float64,1}
end
undrainedaxstresssteps(x) = undrainedaxstresssteps(x,x,x,x)
undrainedaxstresssteps() = undrainedaxstresssteps([0])

struct loadindices
    iloadstart::Array{Int64,1}
    iaxstrain::Array{Int64,1}
    iloadend::Array{Int64,1}
end

struct Pcindices
    iPcstart::Array{Int64,1}
    ilinPcstart::Array{Int64,1}
    ilinPcend::Array{Int64,1}
    iPcend::Array{Int64,1}
end
