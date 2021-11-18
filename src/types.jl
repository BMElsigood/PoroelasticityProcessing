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
