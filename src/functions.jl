
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
    linEν_u(mechStress,mechStrain,iStart,iEnd)
Calculate a combined Young's modulus and Poisson's ratio, (E^u_x)/(1-(ν^u_x))
"""
function linEν_u(mechStress,mechStrain,iStart,iEnd)
    (a,c) = linfit(mechStrain[iStart:iEnd],mechStress[iStart:iEnd])
    Eν = a .*1e-3
    return Eν, c
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
"""
    gassmann(Kd,K0,Kfl,ϕ)
Gassmann fluid substitution to return undrained bulk modulus from drained.
Kd - drained bulk modulus
K0 - bulk modulus of solid
Kfl - fluid modulus
ϕ - porosity
"""
function gassmann(Kd,K0,Kfl,ϕ)
    denom = ϕ / Kfl + (1-ϕ)/K0 -Kd/K0^2
    Ku = Kd + (1 - Kd/K0)^2 / denom
    return Ku
end
"""
    gassmannCijkl(Cd,K0,Kfl,ϕ)
Gassmann fluid substitution to get undrained compliance matrix from drained:
* Cd::Array{Float64,2} - drained compliance matrix
* K0 - bulk modulus of solid
* Kfl - fluid modulus
* ϕ - porosity

### Returns
Cu::Array{Float64,2}
"""
function gassmannCijkl(Cd::Array{Float64,2},K0,Kfl,ϕ)
    Cu = copy(Cd)

    denom = K0 / Kfl * ϕ * (K0 - Kfl) + K0 - (2Cd[1,1] + 2Cd[1,2] + 4Cd[1,3] + Cd[3,3])/9

    Cu11 = Cd[1,1] + (K0 - (Cd[1,1]+Cd[1,2]+Cd[1,3])/3)^2 / denom
    Cu12 = Cd[1,2] + (K0 - (Cd[1,1]+Cd[1,2]+Cd[1,3])/3)^2 / denom
    Cu13 = Cd[1,3] + (K0 - (Cd[1,1]+Cd[1,2]+Cd[1,3])/3)^2 / denom
    Cu33 = Cd[3,3] + (K0 - (Cd[3,3]+2Cd[1,3])/3)^2 / denom

    Cu[1,1] = Cu11
    Cu[1,2] = Cu12
    Cu[1,3] = Cu13
    Cu[2,1] = Cu12
    Cu[2,2] = Cu11
    Cu[2,3] = Cu13
    Cu[3,1] = Cu13
    Cu[2,3] = Cu13
    Cu[3,3] = Cu33
    #Cu44 = Cd44, Cu66 = Cd66

    return Cu
end
"""
    axstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)

Calculates a linear fit of the undrained step in axial stress between indices: istart and iend where istart is the onset of friction
### Returns
DataFrame of length(istart)
columns: |meanstress|Bz|Ez|νz|
"""
function axstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)
    N = length(istart)
    array = zeros(N,4)
    for i in 1:N
        array[i,1] = mean(mechdata.stress[istart[i]:iend[i]])
        array[i,2] = linBz(mechdata.stress,mechdata.pp,istart[i],iend[i])
        array[i,3] = linEz(mechdata.stress,mechdata.εz,istart[i],iend[i])
        array[i,4] = linνz(mechdata.εz,mechdata.εx,istart[i],iend[i])
    end
    colnames = ["meanstress","Bz","Ez","νz"]
    results = DataFrame(array,colnames)
    return results
end
"""
    radstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)

Calculates a linear fit of the undrained step in radial stress (from increasing Pc pump) between indices istart and iend which is the linear part of the increase
### Returns
DataFrame of length(istart)
columns: |meanstress|Bx|Eνz|
"""
function radstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)
    N = size(istart)
    array = zeros(N,3)
    for i in 1:N
        array[i,1] = mean(mechdata.stress[istart[i]:iend[i]])
        array[i,2] = linBx(mechdata.stress,mechdata.pp,istart[i],iend[i])
        array[i,3] = linEνx(mechdata.stress,mechdata.εz,istart[i],iend[i])
    end
    colnames = ["meanstress","Bx","Eνx"]
    results = DataFrame(array,colnames)
    return results
end
