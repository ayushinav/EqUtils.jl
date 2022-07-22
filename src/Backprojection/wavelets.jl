
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