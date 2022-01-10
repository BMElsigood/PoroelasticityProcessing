
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
    linBx(mechStress,mechPp,iStart::Int64,iEnd::Int64)

Calculate Skempton's coefficient Bx from a linear fit of radial stress (Pc at load control) and pore pressure between indices iStart and iEnd
"""
function linBx(mechStress,mechPp,iStart::Int64,iEnd::Int64)
    (a,c) = linfit(mechStress[iStart:iEnd],mechPp[iStart:iEnd])
    Bx = 3a/2
    return Bx
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
    Ez(mechStress,mechStrain,iStart,iEnd)
Calculate Young's modulus between 2 points iStart and iEnd. Used to calculate drained Young's modulus for an undrained step.
"""
function Ez(mechStress,mechStrain,iStart,iEnd)
    ΔStress = mechStress[iEnd] -
                    mechStress[iStart]
    ΔStrain = mechStrain[iEnd] -
                    mechStrain[iStart]
    E = ΔStress ./ ΔStrain
    return E .*1e-3 #GPa
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
    νz(mechStrainz,mechStrainx,iStart,iEnd)

Calculate Poisson's ratio between 2 points iStart, iEnd. Used for drained Poisson's ratio following undrained step.
"""
function νz(mechStrainz,mechStrainx,iStart,iEnd)
    dνradial = mechStrainx[iEnd]-mechStrainx[iStart]
    dνaxial = mechStrainz[iEnd]-mechStrainz[iStart]
    νz = -dνradial/dνaxial
    return νz
end

"""
    linEν_u(mechStress,mechStrain,iStart,iEnd)
Calculate radial stress over strain for radial stress step (Pc pump at load),a combined Young's modulus and Poisson's ratio, (E^u_x)/(1-(ν^u_x)) = 1/ (S11 + S22)
"""
function linEν_u(mechStress,mechStrain,iStart,iEnd)
    (a,c) = linfit(mechStrain[iStart:iEnd],mechStress[iStart:iEnd])
    Eν = a .*1e-3
    return Eν
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
    Cu13 = Cd[1,3] + (K0 - (Cd[1,1]+Cd[1,2]+Cd[1,3])/3)*(K0 - (Cd[3,3]+2Cd[1,3])/3) / denom
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
    sayers95(S::Array{Float64,2},E0,ν0)
Uses Sayers and Kachanov (1995) (doi:10.1029/94JB03134) to get normalised crack density tensors, 2 α and 3 β
Input S is 6x6 stiffness matrix
### Returns
array of noramlised α11,α33,β1111,β1133,β3333
"""
function sayers95(S::Array{Float64,2},E0,ν0)
        C0 = stiffmatrixTI(E0,ν0)
        S0 = inv(C0)
        ΔS1111 = S[1,1] - S0[1,1] #eq 16
        ΔS3333 = S[3,3] - S0[3,3] #eq 17
        ΔS1212 = .25 .*(S[6,6] - S0[6,6]) #eq 18
        ΔS2323 = .25 .*(S[4,4] - S0[4,4]) #eq 19
        ΔS1122 = S[1,2] - S0[1,2] #eq 20
        ΔS2233 = S[2,3] - S0[2,3] #eq 21

        ΔS = [ΔS1111;ΔS3333;ΔS1212;ΔS2323;ΔS1122;ΔS2233]

        A = [1 0 1 0 0; #eq 16
                0 1 0 0 1; #eq 17
                .5 0 .5 0 0; #eq 18
                .25 .25 0 1 0; #eq 19
                0 0 1/3 0 0; #eq 20
                0 0 0 1 0] #eq 21

        αβ = A\ΔS
        γ = αβ .* 3*E0*(2-ν0)/(32*(1-ν0^2))

        return γ
end
"""
    Bx(S,ν0,E0,ϕ,βf)
"""
Bx(S,ν0,E0,ϕ,βf) = 3(S[1,1]+S[1,2]+S[1,3] - (1 - 2ν0)/E0)/
                    (2S[1,1]+S[3,3]+2S[1,2]+4S[1,3] + ϕ*βf -
                        (1+ϕ)*3*(1-2ν0)/E0)
"""
    Bz(S,ν0,E0,ϕ,βf)
"""
Bz(S,ν0,E0,ϕ,βf) = 3(2S[1,3]+S[3,3] - (1 - 2ν0)/E0)/
                    (2S[1,1]+S[3,3]+2S[1,2]+4S[1,3] + ϕ*βf -
                        (1+ϕ)*3*(1-2ν0)/E0)
"""
    B(S,ν0,E0,ϕ,βf)
"""
B(S,ν0,E0,ϕ,βf) = 2/3 .*Bx(S,ν0,E0,ϕ,βf) .+1/3 .*Bz(S,ν0,E0,ϕ,βf)
"""
    BxBzWong(a1,a3,ϕ,βf,E0,ν0)
    Uses Wong (2017) (https://doi.org/10.1002/2017JB014315)

    ### Requires
    * a1,a3 - (normalised) 1st order crack density tensor from Sayers and Kachanov 1995. i.e output [1] and [2] from function Sayers95.
    * ϕ - sample porosity
    * βf - fluid compressibility
    * E0,ν0 - solid material parameters, e.g. from VRH

    ### Returns
    Bx, Bz
"""
function BxBzWong(a1,a3,ϕ,βf,E0,ν0)
    #we want a1 and a3 not normalised
    a1 = a1 ./ (3*E0*(2-ν0)/(32*(1-ν0^2)))
    a3 = a3 ./ (3*E0*(2-ν0)/(32*(1-ν0^2)))
    Bx = a1 /((2*a1 + a3) + ϕ*(βf - 3*(1-2ν0)/E0))
    Bz = a3 /((2*a1 + a3) + ϕ*(βf - 3*(1-2ν0)/E0))
    return Bx,Bz
end
"""
    density(mass,length,diameter)
"""
function density(mass,length,diameter)
    volume = length * π*(diameter/2)^2
    ρ = mass / volume #kg/m3

    return ρ
end
"""
    axstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},idrain::Array{Int64,1},mechdata::keymechparams)

Calculates a linear fit of the undrained step in axial stress between indices: istart and iend where istart is the onset of friction. And drained Young's modulus and Poisson's ratio between 2 points: onset of friction and fully drained.
### Returns
DataFrame of length(istart)
columns: |meanstress|Bz|Ezu|νzu|Ezd|νzd|
"""
function axstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},idrain::Array{Int64,1},mechdata::keymechparams)
    N = length(istart)
    array = zeros(N,6)
    for i in 1:N
        array[i,1] = mean(mechdata.stress[istart[i]:iend[i]])
        array[i,2] = linBz(mechdata.stress,mechdata.pp,istart[i],iend[i])
        array[i,3] = linEz(mechdata.stress,mechdata.εz,istart[i],iend[i])
        array[i,4] = linνz(mechdata.εz,mechdata.εx,istart[i],iend[i])
        array[i,5] = Ez(mechdata.stress,mechdata.εz,istart[i],idrain[i])
        array[i,6] = νz(mechdata.εz,mechdata.εx,istart[i],idrain[i])
    end
    colnames = ["meanstress","Bz","Ezu","νzu","Ezd","νzd"]
    results = DataFrame(array,colnames)
    return results
end
"""
    radstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)

Calculates a linear fit of the undrained step in radial stress (from increasing Pc pump) between indices istart and iend which is the linear part of the increase
### Returns
DataFrame of length(istart)
columns: |meanstress|Bx|Eνx|
where Eνx = (E^u_x)/(1-(ν^u_x)) = 1/ (S11 + S22) (rad stress / rad strain)
"""
function radstresssteps(istart::Array{Int64,1},iend::Array{Int64,1},mechdata::keymechparams)
    N = length(istart)
    array = zeros(N,3)
    for i in 1:N
        array[i,1] = mean(mechdata.stress[istart[i]:iend[i]])
        array[i,2] = linBx(mechdata.Pc,mechdata.pp,istart[i],iend[i])
        array[i,3] = linEν_u(mechdata.Pc,mechdata.εx,istart[i],iend[i])
    end
    colnames = ["meanstress","Bx","Eνx"]
    results = DataFrame(array,colnames)
    return results
end
"""
    axstressstepsdiff(y,istart::Array{Int64,1},ipmax::Array{Int64,1},ipexpundrain::Array{Int64,1},mechdata::keymechparams;
                            Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                            η=0.9096e-3,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15)

Diffusion fit to get permeability, storage capacity
### Returns
DataFrame of length(istart)
columns: |meanstress|perm|stor|Bz|
"""
function axstressstepsdiff(y,istart::Array{Int64,1},ipmax::Array{Int64,1},ipexpundrain::Array{Int64,1},mechdata::keymechparams;
                            Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                            η=0.9096e-3,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15)
    N = length(istart)
    array = zeros(N,4)
    for i in 1:N
        diff = pressurerampsolution(mechdata,y,istart[i],ipmax[i],ipexpundrain[i],Srange=Srange,krange=krange,Brange=Brange,η=η,A=A,L=L,βres=βres)
        array[i,1] = mean(mechdata.stress[istart[i]:ipexpundrain[i]])
        array[i,2] = diff[1]
        array[i,3] = diff[2]
        array[i,4] = diff[3]
    end
    colnames = ["meanstress","perm","stor","Bz"]
    results = DataFrame(array,colnames)
    return results
end

"""
    radstressstepsdiff(y,istart::Array{Int64,1},ipmax::Array{Int64,1},ipexpundrain::Array{Int64,1},mechdata::keymechparams;
                            Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                            η=0.9096e-3,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15)

Diffusion fit to get permeability, storage capacity
### Returns
DataFrame of length(istart)
columns: |meanstress|perm|stor|Bx|
"""
function radstressstepsdiff(y,istart::Array{Int64,1},ipmax::Array{Int64,1},ipexpundrain::Array{Int64,1},mechdata::keymechparams;
                            Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                            η=0.9096e-3,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15)
    N = length(istart)
    array = zeros(N,4)
    for i in 1:N
        diff = pressurerampsolution(mechdata,y,istart[i],ipmax[i],ipexpundrain[i],axial=0,Srange=Srange,krange=krange,Brange=Brange,η=η,A=A,L=L,βres=βres)
        array[i,1] = mean(mechdata.stress[istart[i]:ipexpundrain[i]])
        array[i,2] = diff[1]
        array[i,3] = diff[2]
        array[i,4] = diff[3]
    end
    colnames = ["meanstress","perm","stor","Bx"]
    results = DataFrame(array,colnames)
    return results
end

"""
    lininterp(x::AbstractVector,y::AbstractVector,xi::AbstractVector)
Taken from NBrantut BaseTools
"""
function lininterp(x::AbstractVector,y::AbstractVector,xi::AbstractVector)
    yi = zeros(eltype(y), length(xi))

    # first check if x is sorted and there are no repeated values
    if ~allunique(x)
        error("Values in x should be distinct")
    end

    if ~issorted(x)
        I = sortperm(x)
        x = x[I] #this apparently creates a copy of x, so output not modifed
        y = y[I]
    end

    # sort xi
    if ~issorted(xi)
        J = sortperm(xi)
        xi = xi[J]
    else
        J = 1:length(xi)
    end


    N = length(x)

    for (k,xk) in enumerate(xi)
        m = searchsortedlast(x, xk)
        n = min(max(1, m), N-1)
        yi[J[k]] = y[n] + (xk - x[n])*(y[n+1] - y[n])/(x[n+1]- x[n])
    end

    return yi

end
"""
    lininterp(x::AbstractVector,y::AbstractVector,xi::Number)
Taken from NBrantut BaseTools
"""
function lininterp(x::AbstractVector,y::AbstractVector,xi::Number)

    # first check if x is sorted and there are no repeated values
    if ~allunique(x)
        error("Values in x should be distinct")
    end

    if ~issorted(x)
        I = sortperm(x)
        x = x[I] #this apparently creates a copy of x, so output not modifed
        y = y[I]
    end

    N = length(x)

    m = searchsortedlast(x, xi)
    n = min(max(1, m), N-1)
    yi = y[n] + (xi - x[n])*(y[n+1] - y[n])/(x[n+1]- x[n])

    return yi

end
