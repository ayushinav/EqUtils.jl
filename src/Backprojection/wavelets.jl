
"""Returns the value of ricker wavelet of frequncy `fp` Hz at the instance of `t` sec"""
function rick(fp, t)
    ωp= 2π*fp;
    r= (1 - 0.5*ωp*ωp*t*t) *exp(-0.25*ωp*ωp*t*t);
    return r
end


"""Shifts the wavelet defined by the function `fn` over `tgrid` by a delay of `δt` sec

Positive shift means delay"""
function shifted_wavelet(fn::Function, tgrid, δt::Real)
    return fn.(tgrid.-δt);
end


"""Shifts the wavelet defined by an array `wv` by a delay of `δt/step(tgrid)` where `δt` is in sec"""
function shifted_wavelet(wv::AbstractArray{Float}, tgrid, del_t::Real)
    sh= Int(round(δt/(tgrid[2]- tgrid[1])));
    sh_wv= zeros(size(wv));
    if sh>=0
        for i in 1:length(wv)-sh sh_wv[i+sh]= wv[i] end
    else
        for i in 1:length(wv)-sh sh_wv[i]= wv[i+sh] end
    end
    return sh_wv;
end

"Returns the wavelet contained in the csv data contained in a path (temporary function, needs to be  )"
function path_wav(t)
    if abs(t)>200 return 0
    else
        path_data= readdlm("data_synth_muted.csv", ',');
        path_data= path_data[2:end, 1];
        return path_data[Int(round(2t+ 401))]
    end
end

"Delta function"
delta(t)= if t==0. 1. else 0. end;