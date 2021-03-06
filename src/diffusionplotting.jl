using NumericalIntegration,Roots, SpecialFunctions, QuadGK, LaTeXStrings,PyPlot

"""
    plotcheckfit(mechdata::keymechparams,istart,iend,idrain,Srange,krange,Brange,L;axial = 1,A = π*(40.0e-3 /2)^2,len = 100.0e-3,η=0.9096e-3,βres = 9e-15)

plots "raw" pore pressure verses time, against pp from inversion
"""
function plotcheckfit(mechdata::keymechparams,istart,iend,idrain,Srange,krange,Brange,L;axial = 1,A = π*(40.0e-3 /2)^2,len = 100.0e-3,η=0.9096e-3,βres = 9e-15,nintp = 30)
    I = argmax(exp.(-L))
    Bbest = Brange[I[3]]
    Sbest = Srange[I[2]]
    kbest = krange[I[1]]

    taubest = len^2 * Sbest * η / kbest
    lbest = A*len*Sbest/βres
    #interpolate pp on time for 2 grids
    t0 = mechdata.time[istart:idrain] .- mechdata.time[istart]
    t1 = mechdata.time[istart:iend] .- mechdata.time[istart]
    t2 = mechdata.time[iend+1:idrain] .- mechdata.time[istart]

    t = vcat(t1,
            collect(range(t2[1], stop=t2[end],length=min(nintp,idrain-iend))))

    Pf = mechdata.pp[istart:idrain] .- mechdata.pp[istart]
    Pfnint = lininterp(t0, Pf, t)

    if axial == 0
        ramp = linfit(mechdata.time[istart:iend],mechdata.Pc[istart:iend])[1]
        Bi = round(3/2*Bbest,sigdigits=3)
        txt = L"B_x = "*"$Bi"
    elseif axial == 1
        ramp = linfit(mechdata.time[istart:iend],mechdata.stress[istart:iend])[1]
        Bi = round(3*Bbest,sigdigits=3)
        txt = L"B_z = "*"$Bi"*
            string("\nequilibrate =",round(ramp*Bbest*(mechdata.time[iend]- mechdata.time[istart])*lbest/(lbest+2),sigdigits=3))
    else
        ramp = linfit(mechdata.time[istart:iend],mechdata.Pc[istart:iend])[1]
        Bi = round(Bbest,sigdigits=3)
        txt = L"B = "*"$Bi"
    end

    pcalc = p_ramp(t,taubest, 0.0, lbest,Bbest,ramp,mechdata.time[iend]- mechdata.time[istart],rootsf(lbest, 10)[2:end],len)

    figure()
    xlabel("Time,s")
    ylabel("Pp change, MPa")
    plot(t,pcalc,"r.")
    plot(t0,Pf,"k.")
    axvline(mechdata.time[iend] - mechdata.time[istart])
    annotate(txt,[0.5,0.95],xycoords="axes fraction")

end

function differrors(Srange,krange,Brange,L;nxticks=3,axial=1)
    Sstart = log10(Srange[1])
    Send = log10(Srange[end])
    Slen = length(Srange)

    kstart = log10(krange[1])
    kend = log10(krange[end])
    klen = length(krange)

    Bstart = Brange[1]
    Bend = Brange[end]
    Blen = length(Brange)

    pdf_B = trapz(log.(krange),[trapz(log.(Srange),exp.(-L[i,:,:])) for i in 1:length(krange)])
    pdf_S = trapz(log.(krange),[trapz(Brange,permutedims(exp.(-L[i,:,:]))) for i in 1:length(krange)])
    pdf_k = trapz(Brange,[trapz(log.(Srange),permutedims(exp.(-L[:,:,i]))) for i in 1:length(Brange)])

    pdf_B = pdf_B ./ sum(pdf_B)
    pdf_S = pdf_S ./ sum(pdf_S)
    pdf_k = pdf_k ./ sum(pdf_k)

    Bmult = 3
    if axial == 0
        Bmult = 3/2
    end
    #=
    #95% confidence (2-tailed)
    B95l = lininterp(cumsum(pdf_B),Brange .*Bmult,0.025)
    B95h = lininterp(1 .-cumsum(pdf_B),Brange .*Bmult,0.025)

    #90% confidence (2-tailed)
    B90l = lininterp(cumsum(pdf_B),Brange .*Bmult,0.05)
    B90h = lininterp(1 .-cumsum(pdf_B),Brange .*Bmult,0.05)

    #75% confidence (2-tailed)
    B75l = lininterp(cumsum(pdf_B),Brange .*Bmult,0.125)
    B75h = lininterp(1 .-cumsum(pdf_B),Brange .*Bmult,0.125)

    #50% confidence (2-tailed)
    B50l = lininterp(cumsum(pdf_B),Brange .*Bmult,0.25)
    B50h = lininterp(1 .-cumsum(pdf_B),Brange .*Bmult,0.25)
=#
    I = argmax(exp.(-L))

    figure()
    subplot(311)
    plot(Brange .*Bmult,pdf_B,"k")
#=
    axvline(B95l,c="k")
    axvline(B95h,c="k")
    axvline(B90l,c="k",alpha=0.9)
    axvline(B90h,c="k",alpha=0.9)
    axvline(B75l,c="k",alpha=0.75)
    axvline(B75h,c="k",alpha=0.75)
    axvline(B50l,c="k",alpha=0.5)
    axvline(B50h,c="k",alpha=0.5)
=#
    axvline(Bmult .*Brange[I[3]],c="r")
    if axial ==1
        xlabel("Bz")
    else xlabel("Bx")
    end

    subplot(312)
    semilogx(Srange,pdf_S,"k")
    axvline(Srange[I[2]],c="r")
    xlabel("Storage, 1/Pa")
    #xticks(ticks=range(0,Slen-1,length=nxticks),
    #        labels=round.(10 .^range(Sstart,stop=Send,length=nxticks),sigdigits=3))

    subplot(313)
    semilogx(krange,pdf_k,"k")
    axvline(krange[I[1]],c="r")
    xlabel("Permeability, "* L"m^2")
    #xticks(ticks=range(0,klen-1,length=nxticks),
    #        labels=round.(10 .^range(kstart,stop=kend,length=nxticks),sigdigits=3))
    tight_layout()


end
