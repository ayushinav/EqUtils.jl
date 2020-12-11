
function trace_normalize!(dobs)
    norms=[]
    for ir in 1:size(dobs,2)
        d=view(dobs,:,ir)
        k=inv(norm(d,2))
        push!(norms, k)

        rmul!(d,k)
    end
    return norms
end


# return normalized autocorrelation

function normalized_autocorr(data_obs, dt; slags=[size(data_obs,1)-1, size(data_obs,1)-1])
    nt=size(data_obs,1)
    nr=size(data_obs,2)
    trace_normalize!(data_obs)
    p_conv=Conv.P_conv(eltype(data_obs),gsize=[nt,nr], dsize=[nt,nr], ssize=[sum(slags)+1,nr], slags=slags)
    Conv.mod!(p_conv, d=data_obs, g=data_obs, Conv.S());
    sautos=p_conv.s
    tgridxcorr=dt*range(-slags[1], stop=slags[2], step=1)
    return tgridxcorr, sautos
end


# return envelope of normalized auto-correlations
function normalized_Eautocorr(data_obs,dt)
    txgrid, sautos=normalized_autocorr(data_obs,dt)
    return normalized_E(sautos, dt)
end

# return envelope 
function normalized_E(sautos,dt)
    nt=size(sautos,1)
    ntsp=div(nt,2)+1
    nr=size(sautos,2)
    Esautos=zeros(ntsp,size(sautos,2))
    l2s=zeros(size(sautos,2))
    for i in 1:size(sautos,2)
        s=view(sautos,:,i)
        sss=view(Esautos,:,i)
        ss=copy(s)
	copyto!(sss, abs.(hilbert(ss))[ntsp:end])
        l2s[i]=norm(sss,2)
    end
    return Esautos, l2s
end


# return power specturm using input sautos 
function get_samps(sautos, dt)
    fqs=rfftfreq(size(sautos,1), inv(dt))
    periods=inv.(fqs)
    samps1=abs.(rfft(sautos, [1]));
    return fqs, samps1
end

# get source time functions, assuming a minimum phase
function get_minimum_phase_sources(sautos, dt; env_flag=true)

    fqs=rfftfreq(size(sautos,1), inv(dt))
    periods=inv.(fqs)
    samps1=abs.(rfft(sautos, [1]));

    samps1=sqrt.(samps1)

    sources=[]
    for i in 1:size(samps1,2)
	    S=minimum_phase_signal(samps1[:,i])
	    s=real.(ifft(S))
	    if(env_flag)
		senv=abs.(hilbert(s))
	    	push!(sources, senv)
	    else
	    	push!(sources, s)
	    end
    end
    return hcat(sources...)
end



function get_sauto_specbands(sautos, fmax, dt)
    spec=abs.(rfft(sautos,1))
    nt=size(sautos,1)
    return get_specbands(spec, fmax, nt, dt)
end

function get_specbands(spec, fmax, nt, dt)
	period_bands=[inv(fmax)*0.5^(i-1) for i in 1:4] 
	println(inv.(period_bands))
	#period_bands=[1000, 40, 20, 10]
	f=DSP.rfftfreq(nt, inv(dt))
    period_bands=[1e10, inv(fmax), inv(maximum(f))]

    #nfwidth=round(Int,freq_window_width*inv(step(f)))
    a=zeros(size(spec,2),length(period_bands)-1)
    for ib in 1:length(period_bands)-1
        for i in 1:size(spec,2)
            i1=argmin(abs.(f.-inv(period_bands[ib])))
            i2=argmin(abs.(f.-inv(period_bands[ib+1])))
            ss=view(spec,i1:i2,i)
            a[i,ib]=sum(ss)
#println(i1, "\t",i2, size(f), "\t",f[end])
        end
        aa=view(a,:,ib)
        # normalize!(aa)
    end
#    specbands=copy(a)
    return a, period_bands
end


#function plot_specbands(θ, specbands)
#    p=plot((θ),(specbands[:,1]),proj=:polar, c=:red)
#    plot!(p,(θ),(specbands[:,2]),proj=:polar, c=:green)
#    plot!(p,(θ),(specbands[:,3]),proj=:polar, c=:blue)
#    return p
#end



"""
input amplitude spectrum to return the complex spectrum that has minimum phase
"""
function minimum_phase_signal(F)
	n=length(F)
	

	# to prevent log(0)
	damp = 1e-20 * maximum(abs.(F))
	# logarithm 
	X = log.(complex.(abs.(F) .+ damp, 0.0)) 
	# to cepstral domain - IFFT
	ifft!(X)
	# only real part output
	X = complex.(real.(X), 0.0);
	# scaling
	# X = X / complex(real(npow2), 0.0)

	# positive cepstrum x 2
	X[2 : div(n,2) + 1] .*= complex(2.0, 0.0)
	# remove negative quefrencies
	X[div(n,2) + 2 : n] .= complex(0.0, 0.0) 

	# FFT
	fft!(X)

	# exponential
	F = exp.(X)

end

