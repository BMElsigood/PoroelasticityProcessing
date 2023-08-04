using ProgressMeter

"""
	psolvect(keymechdata::,y,istart::Array{Int64,1},iend::Array{Int64,1},idrain::Array{Int64,1},Srange,krange,Brange;axial = 1,A = π*(40.0e-3 /2)^2,L = 100.0e-3,βres = 9e-15,nintp = 30)

solves pore pressure diffusion equation for step included in iarrays and averages 3D likelyhood matrix

### Optional arguments
axial = 1 (default) for axial stress steps, 0 for radial
A = π*(40.0e-3 /2)^2
L = 100.0e-3
βres = 9e-15
nintp = 30 (default) #resample time using all of ramp points and nintp when diffusing

### Returns
krange,Srange,Brange,average(likelyhoodmatrix)
"""
function psolvect(mechdata::keymechparams, y, istart::Array{Int64, 1}, iend::Array{Int64, 1}, idrain::Array{Int64, 1}, Srange, krange, Brange
	; axial = 1, A = π * (40.0e-3 / 2)^2, len = 100.0e-3, βres = 9e-15, nintp = 30)
	#first check equal length of istart,iend,idrain
	length(istart) == length(iend) == length(idrain) || throw(DimensionMismatch("length of istart, iend, idrain not equal"))

	N = length(istart)
	Lvect = []
	for i in 1:N
		perm, stor, B, L = pressurerampsolution(mechdata, y, istart[i], iend[i], idrain[i], axial = axial, A = A, L = len, Srange = Srange, krange = krange, Brange = Brange, βres = βres, nintp = 30)
		push!(Lvect, L)
	end

	return Srange, krange, Brange, sum(Lvect) ./ N
end
"""
pressurerampsolution(mechdata::keymechparams,y,ipstart,ipmax,ipexpundrain;
								axial = 1,pswitch=0,Srange=10 .^range(-12,stop=-9,length=20),krange=10 .^range(-20,stop=-16,length=50),Brange = range(0,stop=1,length=50),
								η=0.9096e-3,A = π*(40.40e-3 /2)^2,L = 100.9e-3,βres = 9e-15,nintp = 30)
### Optional arguments
axial = 1 (default) for axial stress steps, 0 for radial
pswitch = 0 (default), no plotting, 1 plotting
η=0.9096e-3 #Pa s
A = π*(40.0e-3 /2)^2
L = 100.0e-3
βres = 9e-15
nintp = 30 (default) #resample time using all of ramp points and nintp when diffusing
"""
function pressurerampsolution(mechdata::keymechparams, y, ipstart, ipmax, ipexpundrain;
	axial = 1, Srange = 10 .^ range(-12, stop = -9, length = 20), krange = 10 .^ range(-20, stop = -16, length = 50), Brange = range(0, stop = 1, length = 50),
	η = 0.9096e-3, A = π * (40.0e-3 / 2)^2, L = 100.0e-3, βres = 9e-15,
	nintp = 30)

	# range of exploration for grid search
	ℓrange = [A * L * S / βres for S in Srange]
	τℓrange = [L * η * βres / k / A for k in krange]

	# re-zero Pf for each interval
	pstart = mechdata.pp[ipstart]
	pmax = mechdata.pp[ipmax]
	pexpundrain = mechdata.pp[ipexpundrain]

	Pf = mechdata.pp[ipstart:ipexpundrain] .- pstart

	# interpolate on split time grid (1 during loading, 1 after)
	t0 = mechdata.time[ipstart:ipexpundrain] .- mechdata.time[ipstart]
	t1 = mechdata.time[ipstart:ipmax] .- mechdata.time[ipstart]
	t2 = mechdata.time[ipmax+1:ipexpundrain] .- mechdata.time[ipstart]
	t = vcat(t1, collect(range(t2[1], stop = t2[end], length = min(nintp, ipexpundrain - ipmax))))


	Pfnint = lininterp(t0, Pf, t)
	trampstop = mechdata.time[ipmax] - mechdata.time[ipstart]
	#ramp defined by axial stress change or Pc change
	if axial == 0
		ramp = linfit(mechdata.time[ipstart:ipmax], mechdata.Pc[ipstart:ipmax])[1]
	else
		ramp = linfit(mechdata.time[ipstart:ipmax], mechdata.stress[ipstart:ipmax])[1]
	end
	# fit ℓ and τℓ
	(ℓbest, τℓbest, Bbest, pcalc, likelyhood) = invertp2ramp(t, Pfnint, y, ℓrange, τℓrange, Brange, ramp, trampstop, L)

	perm = (L * βres * η) / (A) ./ τℓbest
	stor = ℓbest * βres / (A * L)
	tau = τℓbest * ℓbest
	l = ℓbest
	B = Bbest

	if axial == 1
		B = 3 .* B #Bz
	else
		B = 3 ./ 2 .* B #Bx
	end
	return perm, stor, B, likelyhood
end

trapz(x, y) = integrate(x, y, Trapezoidal())

# function for which we need roots
function _f_(x, l)
	return 2x + l .* tan.(x)
end

#smart search of roots of f:
function rootsf(l, N)
	out = Array{Float64, 1}(undef, N)

	tol = 1e-12
	#find first root
	x0 = find_zero(x -> _f_(x, l), (0.0, 0.5 * pi - tol), atol = tol)
	out[1] = x0

	c = 1
	while c < N
		x0 = find_zero(x -> _f_(x, l), (0.5 * pi + (c - 1) * pi + tol, 0.5 * pi + c * pi - tol), atol = tol)
		c += 1
		out[c] = x0

	end

	return out
end

function invertp2ramp(t, pobs, y::Float64, lrange, taulrange, Brange, ramp, t0, len)
	sigma = [1]
	return invertp2ramp(t, pobs, y, lrange, taulrange, Brange, sigma, ramp, t0, len)
end

function invertp2ramp(t, pobs, y::AbstractVector, lrange, taulrange, Brange, ramp, t0, len)
	sigma = ones(size(y))
	return invertp2ramp(t, pobs, y, lrange, taulrange, Brange, sigma, ramp, t0, len)
end

function invertp2ramp(t, pobs, y, lrange, taulrange, Brange, sigma, ramp, t0, len)

	NB = length(Brange)
	Nl = length(lrange)
	Nt = length(taulrange)
	L = zeros(Float64, Nt, Nl, NB)

	@showprogress for (iB, B) in enumerate(Brange)
		for (il, l) in enumerate(lrange)
			phi = rootsf(l, 10)
			for (it, taul) in enumerate(taulrange)
				tau = taul * l
				pcalc = p_ramp(t, tau, y, l, B, ramp, t0, phi[2:end], len)
				L[it, il, iB] = sum(sum(abs.(pcalc .- pobs), dims = 2) ./ permutedims(sigma))
			end
		end
	end

	I = argmax(exp.(-L))
	lbest = lrange[I[2]]
	taulbest = taulrange[I[1]]
	taubest = taulbest * lbest
	Bbest = Brange[I[3]]

	pcalc = p_ramp(t, taubest, y, lbest, Bbest, ramp, t0, rootsf(lbest, 10)[2:end], len)

	return lbest, taulbest, Bbest, pcalc, L
end

# with ramp of stress A = (Bz/3) dσz/dt or A = B dPc/dt, until time t0 when stress is constant
function p_ramp(t::Float64, tau, y::Float64, l, B, ramp, t0, phi, len)
	A = B * ramp * tau
	if t < t0
		r = A * (6 + l + 12 * l * (2 + l) * t / tau - 12 * (2 + l) * (y / len)^2) / (12 * (2 + l)^2)
		for phik in phi
			r += (A / (phik^2)) * exp(-4 * phik^2 * t / tau) * cos(2 * phik * y / len) / ((2 + l) * cos(phik) - 2 * phik * sin(phik))
		end
	else
		r = A * t0 / tau * (l / (l + 2))
		for phik in phi
			r += (A / (phik^2)) * (exp(-4 * phik^2 * t / tau) - exp(-4 * phik^2 * ((t - t0) / tau))) * cos(2 * phik * y / len) / ((2 + l) * cos(phik) - 2 * phik * sin(phik))
		end
	end
	return r
end

function p_ramp(t::AbstractVector, tau, y::Float64, l, B, ramp, t0, phi, len)
	return map(x -> p_ramp(x, tau, y, l, B, ramp, t0, phi, len), t)
end

function p_ramp(t::AbstractVector, tau, y::AbstractVector, l, B, ramp, t0, phi, len)
	out = zeros(length(t), length(y))
	for (i, yy) in enumerate(y)
		out[:, i] = p_ramp(t, tau, yy, l, B, ramp, t0, phi, len)
	end
	return out
end
