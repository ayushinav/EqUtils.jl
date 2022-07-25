# Point source
"""`source!(data, mgrid, rgrid, v0, mcntr, izloc::T, iyloc::T, ixloc::T, iTloc::T) where {T<:Integer}` </br>
Updates `data` when a point source were placed at `zgrid`[`izloc`], `ygrid`[`iyloc`], `xgrid`[`ixloc`], and firing at `Tgrid`[`iTloc`].

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**
`mcntr`= Center of the source grids, should ideally by [0,0,0]
`wav`= Source wavelet"""
function source!(data, mgrid, rgrid, v0, mcntr, wav, izloc::T, iyloc::T, ixloc::T, iTloc::T) where {T<:Integer}
    fill!(data, 0.);
    zgrid, ygrid, xgrid, Tgrid= mgrid;
    θ, ϕ= rgrid;
    nr= length(θ)
    zloc= zgrid[izloc]
    yloc= ygrid[iyloc]
    xloc= xgrid[ixloc]
    sloc= [zloc, yloc, xloc]
    dt_fw= Tgrid[iTloc].+ [[sin(θ[i]), sin(ϕ[i])*cos(θ[i]), cos(ϕ[i])*cos(θ[i])]⋅ (sloc.-mcntr)  for i in 1:nr]./v0;
    for ir in 1:nr
        broadcast!(+, view(data, :, ir), view(data, :, ir), gpu(shifted_wavelet(wav, tgrid, dt_fw[ir])));
    end
end

# Extended source
"""`source!(data, mgrid, rgrid, v0, mcntr, izloc::AbstractVector{T}, iyloc::AbstractVector{T}, ixloc::AbstractVector{T}, iTloc::Integer) where {T<:Integer}` </br>
Updates `data` when an extended source were placed at `zgrid`[`izloc`], `ygrid`[`iyloc`], `xgrid`[`ixloc`], all firing simultaneously at `Tgrid`[`iTloc`].

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**
`mcntr`= Center of the source grids, should ideally by [0,0,0]
`wav`= Source wavelet"""
function source!(data, mgrid, rgrid, v0, mcntr, wav, izloc::AbstractVector{T}, iyloc::AbstractVector{T}, ixloc::AbstractVector{T}, iTloc::Integer) where {T<:Integer}
    fill!(data, 0.);
    zgrid, ygrid, xgrid, Tgrid= mgrid;
    θ, ϕ= rgrid;
    nr= length(θ)
    zloc= zgrid[izloc]
    yloc= ygrid[iyloc]
    xloc= xgrid[ixloc]
    slocs= zip(zloc, yloc, xloc)
    ns= length(slocs)
    dt_fw= Tgrid[iTloc].+ [[sin(θ[i]), sin(ϕ[i])*cos(θ[i]), cos(ϕ[i])*cos(θ[i])]⋅ (sloc.-mcntr)  for i in 1:nr, sloc in slocs]./v0;
    for ir in 1:nr, is in 1:ns
        broadcast!(+, view(data, :, ir), view(data, :, ir), gpu(shifted_wavelet(wav, tgrid, dt_fw[ir, is])));
    end
end

# Propagating source
"""`source!(data, mgrid, rgrid, v0, mcntr, izloc::AbstractVector{T}, iyloc::AbstractVector{T}, ixloc::AbstractVector{T}, iTloc::AbstractVector{T}) where {T<:Integer}` </br>
Updates `data` when a propagating source were placed at `zgrid`[`izloc`], `ygrid`[`iyloc`], `xgrid`[`ixloc`], firing at times given by `Tgrid`[`iTloc`].

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**
`mcntr`= Center of the source grids, should ideally by [0,0,0]
`wav`= Source wavelet"""
function source!(data, mgrid, rgrid, v0, mcntr, wav, izloc::AbstractVector{T}, iyloc::AbstractVector{T}, ixloc::AbstractVector{T}, iTloc::AbstractVector{T}) where {T<:Integer}
    fill!(data, 0.);
    zgrid, ygrid, xgrid, Tgrid= mgrid;
    θ, ϕ= rgrid;
    nr= length(θ)
    zloc= zgrid[izloc]
    yloc= ygrid[iyloc]
    xloc= xgrid[ixloc]
    delays= Tgrid[iTloc]
    slocs= zip(zloc, yloc, xloc)
    dt_fw=  [[sin(θ[i]), sin(ϕ[i])*cos(θ[i]), cos(ϕ[i])*cos(θ[i])]⋅ (sloc.-mcntr)  for i in 1:nr, sloc in slocs]./v0;
    for ir in 1:nr, is in 1:ns
        broadcast!(+, view(data, :, ir), view(data, :, ir), gpu(shifted_wavelet(wav, tgrid, dt_fw[ir, is]+ delays[is])));
    end
end