
function get_unilateral_ASTF(θ, θp, cr, c, tgrid, tmax,)
	gamma=(1.0 .- (cr*inv(c) .* cos.(deg2rad.(θ.-θp))))
	cgamma=cr .* inv.(gamma)
	cgamma=inv.(cgamma) .* minimum(cgamma)
	ntmax=argmin(abs.(tgrid.-tmax))
	nt=length(tgrid)
	s=zeros(length(tgrid), length(θ))
	durations=zeros(length(θ))
	for (i,th) in enumerate(θ)
		ss=view(s,:,i)
		# apparent durantion
		nta=(round(Int,ntmax*cgamma[i]))
		durations[i]=(nta-1)*step(tgrid)
		sa=DSP.Windows.lanczos(nta)
		@assert nta <= nt
		for j in 1:nta
			ss[j]=sa[j]
		end
	end
	return s, durations
end

