# Pressure normalised by DeltaP
# Length normalised by L
# Time normalised by tau = L^2/alpha
# Nondimensional parameter l = A L beta/C_res (A: sample area, beta: sample storage, C_res: reservoir storage)
# Flux normalised by A L beta DeltaP/tau


"""
pressurerampsolution(mecht,mechpp,mechstress,y,ipstart,ipmax,ipexpundrain;
                                axial = 1,pswitch=0,Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                                η=0.9096e-3,A = π*(40.40e-3 /2)^2,L = 100.9e-3,βres = 9e-15)
### Optional arguments
axial = 1 (default) for axial stress steps, 0 for radial
pswitch = 0 (default), no plotting, 1 plotting
η=0.9096e-3 #Pa s
A = π*(40.0e-3 /2)^2
L = 100.0e-3
βres = 9e-15
"""
function pressurerampsolution(mechdata::keymechparams,y,ipstart,ipmax,ipexpundrain;
                                axial = 1,Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
                                η=0.9096e-3,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15)

    l(A,L,βres,S) = A*L*S/βres
    #tau = L^2 * S * η /k
    taul(A,L,βres,η,k) = L * η * βres /k /A

    # range of exploration for grid search
    ℓrange = l.(A,L,βres,Srange)
    τℓrange = taul.(A,L,βres,η,krange)
    #σp = 1.0*ones(size(Pf,2))

    # normalise Pf for each interval
    pstart = mechdata.time[ipstart]
    pmax = mechdata.pp[ipmax]
    pexpundrain = mechdata.pp[ipexpundrain]

    Pf = mechdata.pp[ipstart:ipexpundrain] .- pstart
    Pf1 = mechdata.pp[ipstart:ipmax+50] .- pstart

    # interpolate on regular time grid
    # pressure data
    t0 = mechdata.time[ipstart:ipexpundrain] .-mechdata.time[ipstart]
    t = collect(range(t0[1], stop=t0[end],length=30))

    Pfnint = lininterp(t0, Pf, t)
    trampstop = mechdata.time[ipmax]-mecht[ipstart]
    if axial == 1
        ramp = linfit(mechdata.time[ipstart:ipmax],mechdata.stress[ipstart:ipmax])[1]
    else ramp = linfit(mechdata.time[ipstart:ipmax],mechdata.Pc[ipstart:ipmax])[1]
    end
    # fit ℓ and τℓ
    (ℓbest, τℓbest,Bbest, pcalc, likelyhood) = invertp2ramp(t, Pfnint, y, ℓrange, τℓrange,Brange,ramp,trampstop)

    perm = (L*βres*η)/(A)./τℓbest
    stor = ℓbest*βres/(A*L)
    tau = τℓbest*ℓbest
    l = ℓbest
    B = Bbest
    phi = rootsf(l,40)
    pcalc = p_ramp(t0,tau, y, l,B,ramp,trampstop, phi[2:end])

    if axial == 1
        B = 3 .*B #Bz
    else B = 3 ./2 .*B #Bx
    end
    return perm,stor,B,likelyhood
end

trapz(x, y) = integrate(x,y,Trapezoidal())

# function for which we need roots
function _f_(x,l)
    return 2x + l.*tan.(x)
end

#smart search of roots of f:
function rootsf(l, N)
    out = Array{Float64,1}(undef,N)

    tol = 1e-12
    #find first root
    x0 = find_zero(x->_f_(x,l), (0.0, 0.5*pi-tol), atol=tol)
    out[1] = x0

    c=1
    while c<N
        x0 = find_zero(x->_f_(x,l), (0.5*pi + (c-1)*pi +tol, 0.5*pi + c*pi -tol), atol=tol)
        c+=1
        out[c] = x0

    end

    return out
end

function invertp2ramp(t, pobs, y, lrange, taulrange,Brange,ramp,t0)
    sigma = ones(size(y))
    return invertp2ramp(t, pobs, y, lrange, taulrange,Brange, sigma,ramp,t0)
end

function invertp2ramp(t, pobs, y, lrange, taulrange,Brange, sigma,ramp,t0)

    NB = length(Brange)
    Nl = length(lrange)
    Nt = length(taulrange)
    L = zeros(Float64, Nt, Nl, NB)

    for (iB, B) in enumerate(Brange)
        println("loop: iB = $iB")
        for (il, l) in enumerate(lrange)
            phi = rootsf(l,40)
            for (it, taul) in enumerate(taulrange)
                tau = taul*l
                pcalc = p_ramp(t,tau,y,l,B,ramp,t0,phi[2:end])
                L[it,il,iB] = sum(sum(abs.(pcalc .- pobs),dims=2)./permutedims(sigma))
            end
        end
    end

    println("size(L) = ",size(L))
    I = argmax(exp.(-L))
    lbest = lrange[I[2]]
    taulbest = taulrange[I[1]]
    taubest = taulbest*lbest
    Bbest = Brange[I[3]]

    pcalc = p_ramp(t,taubest, y, lbest,Bbest,ramp,t0,rootsf(lbest, 40)[2:end])
    println("size pcalc = ",size(pcalc))
    #=
    #compute marginal pdfs
    # here we use log of x,y values because they are jeffreys parameters!!
    pdf_taul = trapz(log.(Brange),trapz(log.(lrange), permutedims(exp.(-L))))

    plot(pdf_taul)

    pdf_ell  = trapz(log.(taulrange), exp.(-L))
    # compute "errors" based on marginal pdfs
    i_taul = findall(pdf_taul.>maximum(pdf_taul)/2)
    taulmin = taulrange[i_taul[1]]
    taulmax = taulrange[i_taul[end]]

    i_ell = findall(pdf_ell.>maximum(pdf_ell)/2)
    ellmin = lrange[i_ell[1]]
    ellmax = lrange[i_ell[end]]
    =#

    return lbest, taulbest,Bbest, pcalc, L#, (ellmin, ellmax), (taulmin, taulmax)
end

# with ramp of stress A = (Bz/3) dσz/dt or A = B dPc/dt, until time t0 when stress is constant
function p_ramp(t::Float64,tau,y::Float64,l,B,ramp,t0, phi)
    A = B*ramp
    if t < t0
        r = A*t*(1- 2/(l+2))
        for phik in phi
            r += (tau*A/(4*phik^2))*(1-exp(-4*phik^2 *t /tau))*(2.0*cos(2*phik*y/L)*sin(phik))/(phik - cos(phik)sin(phik))
        end
    else
        r = A*t0*(1-2/(l+2))
        for phik in phi
            r += (tau*A/(4*phik^2))*(exp(-4*phik^2 *(t-t0)/tau)-exp(-4*phik^2 *t/tau))*(2.0*cos(2*phik*y/L)*sin(phik))/(phik - cos(phik)sin(phik))
        end
    end
    return r
end

function p_ramp(t::AbstractVector,tau,y::Float64,l,B,ramp,t0, phi)
    return map(x->p_ramp(x,tau, y, l,B,ramp,t0, phi), t)
end

function p_ramp(t::AbstractVector,tau, y::AbstractVector,l,B,ramp,t0, phi)
    out = zeros(length(t), length(y))
    for (i,yy) in enumerate(y)
        out[:,i] = p_ramp(t,tau, yy,l,B,ramp,t0, phi)
    end
    return out
end
