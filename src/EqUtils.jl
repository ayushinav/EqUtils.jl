module EqUtils

using StatsBase
using SAC
using CSV
using Signals
using FocusedBlindDecon
using Dates
using LinearAlgebra
using DataFrames
using MLDataPattern
using DSP
using FFTW
using Distributions
using Random
using StatsBase


include("circ_mean.jl")
include("uniquetol.jl")


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

	i=findall((gcarcs .< g1) | (gcarcs .> g2))

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

function cut_noisy_records!(Tg, filenames, tnoise, threshold) 
	tgrid=range(Tg[1].b, stop=Tg[1].e, length=Tg[1].npts)
	dobs=hcat([(Tg[i].t) for i in 1:length(Tg)]...)
	# energy in first tnoise seconds
	noise=[norm(dobs[1:round(Int,tnoise*inv(step(tgrid))),i]) for i in 1:size(dobs,2)]; 

	icut=[true for i in 1:length(Tg)] # all are ON
	i=findall(noise .> threshold)

	deleteat!(Tg, i)
	deleteat!(filenames, i)

	println(string("removed ",length(i)," noisy receivers\n"))
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

	println("fwewec")
	d1=complex.(gcarcs.*cos.(deg2rad.(azs)),  gcarcs.*sin.(deg2rad.(azs)))


	ii=fill(false, length(Tg))
	for az in azgrid
		for gcarc in gcarcgrid
			d2=complex(gcarc*cos(deg2rad(az)),  gcarc*sin(deg2rad(az)))
			dist=abs.(d1.-d2)
			iii=findall(dist.==minimum(dist))
			for iiii in iii
				ii[iiii]=true
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

	println("Number of Records Removed:\t",count(.!ii) )
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
