"""
    taken from Brantut BaseTools.jl
    linfit(x,y)
"""
function linfit(x, y)
    if length(x)!=length(y)
        error("x and y should be of the same size")
    end

    A = [x  ones(size(x,1))]

    coefs = A\y

    return coefs[1], coefs[2]
end

"""
    linBz(mechStress,mechPp,iStart::Int64,iEnd::Int64)

Calculate Skempton's axial coefficient Bz from a linear fit of axial stress and pore pressure between indices iStart and iEnd
"""
function linBz(mechStress,mechPp,iStart::Int64,iEnd::Int64)
    (a,c) = linfit(mechStress[iStart:iEnd],mechPp[iStart:iEnd])
    Bz = 3a
    return Bz
end

"""
    linB(mechStress,mechPp,iStart::Int64,iEnd::Int64)

Calculate Skempton's coefficient B from a linear fit of confining pressure and pore pressure between indices iStart and iEnd
"""
function linB(mechStress,mechPp,iStart::Int64,iEnd::Int64)
    (a,c) = linfit(mechStress[iStart:iEnd],mechPp[iStart:iEnd])
    B = a
    return B
end

"""
    linEz(mechStress,mechStrain,iStart,iEnd)

Calculate Young's modulus Ez (GPa) from a linear fit of (axial) stress (MPa) and (axial) strain between iStart and iEnd.
"""
function linEz(mechStress,mechStrain,iStart::Int64,iEnd::Int64)
    (a,c) = linfit(mechStrain[iStart:iEnd],mechStress[iStart:iEnd])
    E = a .*1e-3
    return E
end

"""
    volStrain(axStrain,radStrain)
"""
volStrain(axStrain,radStrain) = @. axStrain + 2*radStrain

"""
    linK(mechPc,mechStrainz,mechStrainx,iStart,iEnd)
Calculate bulk modulus in (GPa) using a linear fit of confining pressure (MPa) and volumetric strain
"""
function linK(mechPc,mechStrainz,mechStrainx,iStart,iEnd)
    εvol = volStrain(mechStrainz[iStart:iEnd],mechStrainx[iStart:iEnd])
    (a,c) = linfit(εvol,mechPc[iStart:iEnd])
    K = a
    return K .*1e-3 #GPa
end
"""
    linνz(mechStrainz,mechStrainx,iStart::Int64,iEnd::Int64)
Calculate Poisson's ratio -εx/εz using linfit between iStart and iEnd
"""
function linνz(mechStrainz,mechStrainx,iStart::Int64,iEnd::Int64)
    (a,c) = linfit(mechStrainz[iStart:iEnd],mechStrainx[iStart:iEnd])
    ν = -a
    return ν
end
"""
    darcyperm(t,p1,p2,vol,area,L,μ,iStart,iEnd)
Calculate constant flow permeability where p1 and p2 are different measured pore pressures, length L apart
"""
function darcyperm(t,p1,p2,vol,area,L,μ,iStart,iEnd)
    (a,c) = linfit(t[iStart:iEnd],vol[iStart:iEnd])#dV/dt
    Q = a
    Δp = abs(mean(p1[iStart:iEnd])-mean(p2[iStart:iEnd]))
    k = Q*μ*L / (area * Δp)
    return k
end

function axstresssteps(isteps::Array{Int64,2},mechdata::keymechparams)
    N = size(isteps,1)
    array = undrainedaxstresssteps()
    for i in 1:N
        push!(array.meanstress,mean(mechdata.stress[isteps[i,2]:isteps[i,3]]))
        push!(array.Bz,linBz(mechdata.stress,mechdata.pp,isteps[i,2],isteps[i,3]))
        push!(array.Ez,linEz(mechdata.stress,mechdata.εz,isteps[i,2],isteps[i,3]))
        push!(array.νz,linνz(mechdata.εz,mechdata.εx,isteps[i,2],isteps[i,3]))
        #array[i,:]=[meanstress Bz Ez νz]
    end
    return array
end
