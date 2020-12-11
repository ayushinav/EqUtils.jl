module EqUtils


using StatsBase
using Seis
using CSV
#using Signals
using Dates
using LinearAlgebra
using Interpolations
using ProgressMeter
using DataFrames
using Statistics
using MLDataPattern
using DSP
using FFTW
using Distributions
using Random
using StatsBase
using Conv
using TauPy



include("circ_mean.jl")
include("uniquetol.jl")
include("spec_tools.jl")
include("ambiguity.jl")
include("rupture.jl")

"""
cut records just after minimum(P, Pdiff, PKP)
"""
function cut_P(Tg, reference_time::Float64, twidth=500)
	delta=Tg[1].delta
	nt=length(Tg[1].t)

	evdps=getfield.(Tg,:evdp) 
	gcarcs=getfield.(Tg,:gcarc) 
	 # taking average in case of triplications

	 ptime(evdp, x)=minimum(getfield.(travel_time(evdp, x, ["P", "Pdiff", "PKP"]),:time))  

	 ptimes=[ ptime(evdps[i], gcarcs[i]) for i in 1:length(evdps)];   

	 tcuts=reference_time .+ ptimes .+ twidth 
	 itcuts=round.(Int, tcuts .* inv(delta))

	 itcutmin=minimum(itcuts)
	 itcutmax=maximum(itcuts)

	 TgP=deepcopy(Tg)
	 for (i,T) in enumerate(TgP)
		 T.t=zeros(eltype(T.t), itcutmax)
		 for it in 1:itcuts[i]
			 T.t[it]=Tg[i].t[it]
			 T.npts=length(T.t)
		 end
	 end

	 return TgP
end




#=

"""
cut records at mean(travel_time(PPP), travel_time(SS))
"""
function cut_PPP_SS(Tg, reference_time::Float64)
	delta=Tg[1].delta
	nt=length(Tg[1].t)

	evdps=getfield.(Tg,:evdp) 
	gcarcs=getfield.(Tg,:gcarc) 
	 # taking average in case of triplications

	 function stime(evdp, x)
		 ppp_time=Statistics.mean(getfield.(travel_time(evdp, x, ["PPP"]),:time))  
		 if(x<100)
			 return mean([Statistics.mean(getfield.(travel_time(evdp, x, ["S"]),:time)), ppp_time])
		 else
			s_time_100=Statistics.mean(getfield.(travel_time(evdp, 100., ["S"]),:time))    
			ss_time_100=Statistics.mean(getfield.(travel_time(evdp, 100., ["SS"]),:time))    
			return mean([(Statistics.mean(getfield.(travel_time(evdp, x, ["SS"]),:time)) -  (ss_time_100-s_time_100)) , ppp_time])
		end
	end


	 sstimes=[ stime(evdps[i], gcarcs[i]) for i in 1:length(evdps)];   
	 ppptimes=[Statistics.mean(getfield.(travel_time(evdps[i], gcarcs[i], ["PPP"]),:time)) for i in 1:length(evdps)];   

	 tcuts=reference_time .+ [Statistics.mean([sstimes[i], ppptimes[i]]) for i in 1:length(sstimes)]
	 itcuts=round.(Int, tcuts .* inv(delta))

	 itcutmin=minimum(itcuts)
	 itcutmax=maximum(itcuts)

	 TgPPP=deepcopy(Tg)
	 for (i,T) in enumerate(TgPPP)
		 T.t=zeros(eltype(T.t), itcutmax)
		 for it in 1:itcuts[i]
			 T.t[it]=Tg[i].t[it]
		 end
	 end

	 TgSS=deepcopy(Tg)
	 for (i,T) in enumerate(TgSS)
		 T.t=zeros(eltype(T.t), nt-itcutmin+1)
		 for it in itcuts[i]:nt
			 T.t[it-itcutmin+1]=Tg[i].t[it]
		 end
	 end

	 return TgSS




end


=#

function tt(Tg)
	delta=Tg[1][:delta]
	tt=zeros(length(Tg))
	for i in 1:length(Tg)
		k=TauPy.travel_time(Tg[i][:evdp], Tg[i][:gcarc], ["P"])
		if(length(k)==0)
			k=TauPy.travel_time(Tg[i][:evdp], Tg[i][:gcarc], ["Pdiff"])
		end
		if(length(k)≠0)
			tt[i]=getfield(k[1],:time)
		end
	end
	# cut
	itmax=round.(Int,(500+maximum(tt)).*inv(delta)).+1
	for i in 1:length(Tg)
		itcut=round.(Int,(500+tt[i]).*inv(delta)).+1
		Tg[i].t[itcut+1:end] .= 0.0
		a=Tg[i].t[1:itmax]
		Tg[i].t = a
		Tg[i].npts=length(a)
	end
	return tt
end

"""
tbeg is offset starting maximum(nztime)
"""
function synchronize!(Tg, filenames; tbeg=nothing, toffset=nothing)
	npts=mode(getfield.(Tg,:npts))
	nptsvec=getfield.(Tg,:npts)

	# compute the kztime for each record
	nzjday=getfield.(Tg, :nzjday)
	nzhour=getfield.(Tg, :nzhour)
	nzmin=getfield.(Tg, :nzmin)
	nzsec=getfield.(Tg, :nzsec)
	nzmsec=getfield.(Tg, :nzmsec)
	b=getfield.(Tg, :b)
	e=getfield.(Tg, :e)

	# nztime (in seconds) of first sample 
	nztime=nzjday.*(24*60*60).+nzhour.*(60*60).+nzmin.*(60).+nzsec.+nzmsec.*(1e-3)
	if(tbeg===nothing)
		tbeg=maximum(nztime.+b)
	else
		tbeg=maximum(nztime.+b)+tbeg
	end

	if(toffset===nothing)
		tend=minimum(nztime.+e)
	else
		tend=maximum(nztime.+b).+toffset
		@assert all(tend .< minimum(nztime.+e))
	end
		@assert tbeg<tend

	delta=Tg[1][:delta]
	itbeg_rel=round.(Int,(tbeg.-nztime.-b).*inv(delta)).+1
	itend_rel=round.(Int,(tend.-nztime.-b).*inv(delta))
	itoffset=minimum(itend_rel.-itbeg_rel)
	for (i,T) in enumerate(Tg)
		a=copy(T.t)
		T.t=a[itbeg_rel[i] : itbeg_rel[i]+itoffset]
		T.npts=itoffset+1
		T.b=0
		T.e=itoffset*delta
	end
	return nothing
end



"""
* keep only the first mode of length 
* use it to remove some bad records of shorter time duration
"""
function cut_short_records!(Tg, filenames)
	npts=mode(getfield.(Tg,:npts))
	nptsvec=getfield.(Tg,:npts)

	i=findall(nptsvec .≠ npts)

	deleteat!(Tg, i)
	deleteat!(filenames, i)

	return nothing
end


"""
sort by gcarc
"""
function sort_gcarc!(Tg, filenames)
	gcarcs=getfield.(Tg, :gcarc)
	isort=sortperm(gcarcs)
	Tg.=Tg[isort]
	filenames.=filenames[isort]
	return nothing
end

"""
cut records by distance
* only allow between g1 and g2
"""
function cut_gcarc!(Tg, filenames, g1, g2=Inf)
	gcarcs=getfield.(Tg, :gcarc)

	i=findall((gcarcs .< g1) .| (gcarcs .> g2))

	deleteat!(Tg, i)
	deleteat!(filenames, i)

	return nothing
end

function cut_nan!(Tg, filenames)
	i=findall(broadcast(x->any(isnan.(x)), getfield.(Tg,:t)))

	deleteat!(Tg, i)
	deleteat!(filenames, i)

	return nothing
end

"""
compute the energy in first tnoise seconds and then cut traces with 
func == > cut traces with energy more than certain threshild
func == < cut traces with energy less than certain threshild
"""
function cut_weird_records!(Tg, filenames, tnoise, threshold, func) 
	tgrid=range(Tg[1].b, stop=Tg[1].e, length=Tg[1].npts)
	dobs=hcat([(Tg[i].t) for i in 1:length(Tg)]...)
	# energy in first tnoise seconds
	noise=[norm(dobs[1:round(Int,tnoise*inv(step(tgrid))),i]) for i in 1:size(dobs,2)]; 

	icut=[true for i in 1:length(Tg)] # all are ON
	i=findall(func.(noise, threshold))

	deleteat!(Tg, i)
	deleteat!(filenames, i)

	println(string("removed ",length(i)," weird receivers\n"))

	return nothing
end



"""
delete a few records, such that the azimuths are uniformly distributed
"""
function uniformly_sample!(Tg, filenames; δaz=1, δgcarc=3)

	azs=getfield.(Tg, :az) # in degrees
	gcarcs=getfield.(Tg, :gcarc) # in degrees

	azgrid=range(minimum(azs), stop=maximum(azs), step=δaz)
	gcarcgrid=range(minimum(gcarcs), stop=maximum(gcarcs), step=δgcarc)

	println(string(length(azgrid)*length(gcarcgrid),"\t receivers will be choosen roughly"))
	d1=complex.(gcarcs.*cos.(deg2rad.(azs)),  gcarcs.*sin.(deg2rad.(azs)))

	ii=fill(false, length(Tg))
	for az in azgrid
		for gcarc in gcarcgrid
			d2=complex(gcarc*cos(deg2rad(az)),  gcarc*sin(deg2rad(az)))
			dist=abs.(d1.-d2)
			# choose multiple components in sequence
			iii=findall(dist.==minimum(dist)) # findall for multiple components
			# check which rec is off
			irec=findfirst(.!ii[iii])
			# switch it on
			if(!(irec===nothing))
				ii[iii[irec]]=true
			end
		end
	end

	println("Number of Records Removed:\t",count(.!ii) )
	i=findall(.!ii)


	deleteat!(Tg, i)
	deleteat!(filenames, i)

	return nothing
end


"""
First generate a grid, then select nearest elements, finally delete unselected elements.
"""
function uniformly_sample!(v, value, δ)

	grid=range(minimum(value), stop=maximum(value), step=δ)

	ii=fill(false, length(v))
	for az in grid
		dist=abs.(value.-az)
		iii=findall(dist.==minimum(dist))
		for iiii in iii
			ii[iiii]=true
		end
	end

	println("Number of Elements Removed:\t",count(.!ii) )
	i=findall(.!ii)


	deleteat!(v, i)

	return nothing
end

"""
Same as above, but deleteat! applied to two vectors v and w
"""
function uniformly_sample!(v, w, value, δ)
	uniformly_sample!(v, value, δ)
	uniformly_sample!(w, value, δ)
end





end # module
