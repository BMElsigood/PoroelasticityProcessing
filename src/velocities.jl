"""
    svyCompliancet(survey,ρ,(θp,vp)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}},(θsh,vsh)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}};angerr = 0.05,vperr = 0.1,vsherr=0.2,errlogC = 1)

Calculates the compliance matrix from density and the p and sh wavespeeds and angles by inverting thomsen parameters starting with isotropic compliance matrix C0

# Required Inputs
* survey: an survey index for the velocities required
* ρ: density in kg/m3
* (θp,vp)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}}
    - Θp: array containing angle to verticle in radians of p-wave sensors
    - Vp: array where each element is an array of p-wave velocities (m/s) from a survey at corresponding angles Θp
* (θsh,vsh)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}}
    - Θsh: array containing angle to verticle in radians of sh-wave sensors
    - Vsh: array where each element is an array of sh-wave velocities (m/s) from a survey at corresponding angles Θsh

# Optional Arguments
* angerr: default = 0.05, error on log(angles), same for p and sh
* vperr: default = 0.1, error on log(vperr), same for all angles
* vsherr: default = 0.2, error on log(vsherr), same for all angles

## Returns

* C1: 6x6 compliance matrix calculated from inverted thomsen parameters
* cm: covariance matrix from findthomsen()
* anglep: array of phase angles for p-waves from findthomsen()
* anglesh: array of phase angles for sh-waves from findthomsen()

## Useage
C1,cm,anglep,anglesh = svyCompliancet(survey,ρ,(θp,vp),(θsh,vsh))
"""
function svyCompliancet(survey,ρ,(θp,vp)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}},(θsh,vsh)::Tuple{Array{Float64,1},Array{Array{Float64,1},1}}
                            ;angerr = 0.05,vperr = 0.1,vsherr=0.2)

    ρ = ρ*1e-3

    p = VMeasure.(θp,vp[survey].*1e-3,angerr,vperr);
    sh = VMeasure.(θsh,vsh[survey].*1e-3,angerr,vsherr);
    sv = VMeasure[]

    C0 = stiffmatrixTI(vpvsYoungsPoisson(p[1].value,sh[1].value,ρ))
    ϵ, δ, γ = moduli2thomsen(C0)
    vp0,vs0 = velocities0(C0,ρ)

    m1, anglep, anglesv, anglesh, cm = findthomsen(p,sv,sh,(vp0,10),(vs0,10),(ϵ,10),(δ,10),(γ,10))
    #m1, anglep, anglesv, anglesh, cm = findthomsen(p,sv,sh,vp0,vs0,ϵ, δ, γ)

    C1 = thomsen2moduli(m1[1],m1[2],ρ,m1[3],m1[4],m1[5])
    return C1.*1e9,cm,anglep,anglesh
end
"""
    svyCompliancet(survey,ρ,(θp,vp),(θsh,vsh),C0;angerr = 0.05,vperr = 0.1,vsherr=0.2,errlogC = 1)
Calculates the compliance matrix from density and the p and sh wavespeeds and angles by inverting thomsen parameters starting with input compliance matrix C0 (which can be isotropic)
"""
function svyCompliancet(survey,ρ,(θp,vp),(θsh,vsh),C0;angerr = 0.05,vperr = 0.1,vsherr=0.2,errlogC = 1)

    ρ = ρ*1e-3
    C0 = C0.*1e-9

    p = VMeasure.(θp,vp[survey].*1e-3,angerr,vperr);
    sh = VMeasure.(θsh,vsh[survey].*1e-3,angerr,vsherr);
    sv = VMeasure[]

    ϵ, δ, γ = moduli2thomsen(C0)
    vp0,vs0 = velocities0(C0,ρ)

    m1, anglep, anglesv, anglesh, cm = findthomsen(p,sv,sh,(vp0,100),(vs0,100),(ϵ,100),(δ,100),(γ,100))
    #m1, anglep, anglesv, anglesh, cm = findthomsen(p,sv,sh,vp0,vs0,ϵ, δ, γ)

    C1 = thomsen2moduli(m1[1],m1[2],ρ,m1[3],m1[4],m1[5])
    return C1.*1e9,cm,anglep,anglesh
end

function phasetogroupvp(angphase,C,ρ)
    ϵ, δ, γ = moduli2thomsen(C)
    vp0,vs0 = velocities0(C,ρ)
    vp = vphaseP(angphase,C,ρ)
    #vp = vphaseP(angp,vp0,vs0,ϵ,δ)
    dvvp = dvPdthetav(angphase,vp0/vs0,ϵ,δ)
    angpgr = groupangle(angphase,dvvp)
    vpgroup = vgroup(vp, dvvp)

    return vpgroup
end

function phasetogroupvsh(angphase,C,ρ)
    ϵ, δ, γ = moduli2thomsen(C)
    vp0,vs0 = velocities0(C,ρ)
    angshgr = groupangleSH(angphase,γ)
    vshgroup = vgroupSH(angshgr, vs0, γ)

    return vshgroup
end

function surveytocompliancedf(surveyrange,ρ,(θp,vp),(θsh,vsh))
    N = length(surveyrange)
    array = Array{Any}(nothing,N,11)
    for (n,i) in enumerate(surveyrange)
        array[i,1] = i
        array[i,2] = θp
        array[i,3] = Vp
        array[i,4] = θsh
        array[i,5] = Vsh
        if n == 1
            C1,cm,anglep,anglesh = svyCompliancet(survey,ρ,(θp,vp),(θsh,vsh))
        else C1,cm,anglep,anglesh = svyCompliancet(survey,ρ,(θp,vp),(θsh,vsh),C1)
        end
        array[i,6] = C1
        ϵ, δ, γ = moduli2thomsen(C1)
        array[i,7:9] = [ϵ, δ, γ]
        vpgroup = phasetogroupvp.(anglep,C1,ρ)
        array[i,10] = vpgroup
        vshgroup = phasetopgroupvsh.(anglesh,C1,ρ)
        array[i,11] = vshgroup
    end
    colnames = ["Index","θp","Vp","θsh","Vsh","ϵ","δ","γ","vpgroup","vshgroup"]
    results = DataFrame(array,colnames)
    return results
end
#=
not sure this belongs here
"""
    joinVelMech(mechdata,surveydata)
For each survey the mech data is interpolated
p is array of tuples of pairs of paths
Output of dataframe of mechdata and path velocities for each survey
"""
function joinVelMech(mechdata,surveydata,p)

    offset = Dates.datetime2unix(mechdata[:date])-surveydata[1].t

    tMech = mechdata[:t].-mechdata[:t][2] .+offset #zero mech timeline to second point as it is continuous from then

    tZero = surveydata[1].t
    t = [s.t for s in surveydata]

    vPc = lininterp(tMech,mechdata[:Pc],t .- tZero)
    vStress = lininterp(tMech,mechdata[:stress],t .- tZero)
    vpB = lininterp(tMech,mechdata[:pB],t .- tZero)
    vpA = lininterp(tMech,mechdata[:pA],t .- tZero)
    vpp1 = lininterp(tMech,mechdata[:pp1],t .- tZero)
    vpp2 = lininterp(tMech,mechdata[:pp2],t .- tZero)


    vMech = DataFrame()

    vMech.t = t
    vMech.Pc = vPc
    vMech.stress = vStress
    vMech.pB = vpB
    vMech.pA = vpA
    vMech.pp1 = vpp1
    vMech.pp2 = vpp2

    colnames = ["t","Pc","stress","pB","pA","pp1","pp2"]

    for (i,j) in p
        lb1 = string("v$i","_$j")
        push!(colnames,lb1)
        vMech.lb1 = [s.V[i,j] for s in surveydata]

        lb2 = string("v$j","_$i")
        push!(colnames,lb2)
        vMech.lb2 = [s.V[j,i] for s in surveydata]

        rename!(vMech,colnames)
    end
    return vMech
end
=#
