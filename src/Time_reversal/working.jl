using Plots
using SpecialFunctions
using LinearAlgebra
using FFTW
using GPUArrays
using BenchmarkTools
using Optim
using Statistics
using Test
using ProgressMeter
using HDF5
using Printf
using Flux
using ParallelStencil
@init_parallel_stencil(CUDA, Float64, 3)
using StatsBase
using DelimitedFiles
using LineSearches

Float= Float64;
Int= Int32;
xpu= gpu;

include("wavelets.jl");
include("Synthetics.jl");
include("utils_bp.jl");
include("backproject.jl");
include("utils_optim.jl");

v0= 8.000 # km/s Velocity around the source grid
d0= 10000.0; # km

# Time grid
nt = 800
dt= 0.5
fmax= 0.25/dt
tgrid = Float.(range(-200, stop=200, step= dt)); # time samples from -200, 200 seconds
nt= length(tgrid);

# eq= "chl1";
fl= open("tmp.txt", "r")
eq= read(fl, String)
close(fl)

# Different eqs would have different grids, depending on the extent of rupture
dx= 10; dy= 10; 
limx= 500; limy= 500; # km
dT= 0.5;
Tgrid = range(-50, 200, step= dT);


dobs, k= get_symae_data(eq);

cntr, θ, ϕ, travel_times, network, ids, nevids= obspy_data(eq, "/home/abhinav/workspace/Bp_eq_data_processed", "/mnt/data2/isha/unique_pwin33_traces/");
nr= length(θ);

θ, ϕ, travel_times, network, ids, nr, nevids, IDs= filter_nevids(eq, θ, ϕ, travel_times, network, ids, nr, nevids);
dobs= cpu(dobs)[:,IDs] |>xpu;

bins= range(-180, 180, step= 2);#.*π./180;

cntr, θ, ϕ, travel_times, network, nr, IDs= get_unique_receivers_per_bin(cntr, θ, ϕ, travel_times, network, θ, bins, nr, 1);
nevids= nevids[IDs];
dobs= cpu(dobs)[:,IDs] |>xpu;

# cntr, θ, ϕ, travel_times, nevids, nr, IDs= filter_arrays([8], cntr, θ, ϕ, travel_times, nevids, nr)
# dobs= cpu(dobs)[:,IDs] |>xpu;

# cntr, θ, ϕ, travel_times, nevids, nr, IDs= get_unique_receivers(cntr, θ, ϕ, travel_times, nevids, nr)
# dobs= cpu(dobs)[:,IDs] |>xpu;

# rgrid
rx= zeros(nr);
ry= zeros(nr);
rz= zeros(nr);

for i in 1:nr rz[i], ry[i], rx[i]= d0.*to_zyx(θ[i], ϕ[i]) end

zgrid = Float.([0])
ygrid = Float.(range(-limy, stop=limy, step=dy))
xgrid = Float.(range(-limx, stop=limx, step=dx))
mgrid= [zgrid, ygrid, xgrid, Tgrid];
rgrid= [rz, ry, rx];
nz, ny, nx, nT= length.(mgrid)
Tshift= Int.(round.(Tgrid./dt))
shift= get_shifts(mgrid, rgrid, v0, d0);
ricker(t)= delta(t);

@show nz, ny, nx, nT, nt, nr;

wavg= fill(Float(0.), nt, 1, 1, 1, 1, 1) |>xpu; # ricker
shiftg= fill(Int(0), 1, nr, nz, ny, nx, 1)|>xpu; # time shifts
Tshiftg= fill(Int(0), 1,1,1,1,1,nT)|>xpu; #time shits for temporal domain
@allowscalar copyto!(view(wavg, :, 1, 1, 1, 1, 1), path_wav.(tgrid))# ;FFTW.rfft(r));
@allowscalar copyto!(view(shiftg, 1 , :, :, :, :, 1), permutedims(shift, [4, 1, 2, 3]));
@allowscalar copyto!(view(Tshiftg, 1, 1, 1, 1, 1, :), Tshift);
delay= broadcast(+, shiftg, Tshiftg);
(lb, ub)= extrema(tgrid)./dt .- reverse(extrema(cpu(delay)));

exp_tgrid= Float.(dt.*(lb:ub));
wav0= Float.(ricker.(exp_tgrid)) |>xpu;
dg= zeros(Float, nt, nr) |>xpu;
img= zeros(Float, nz, ny, nx, nT) |>xpu;

mvec= zeros(Float, 1, 1, nz, ny, nx, nT) |>xpu;
dvec= zeros(Float, nt, nr, 1, 1, 1, nT) |>xpu;
g2= zeros(Float, nt, nr, nT) |>xpu;
grvec= zeros(Float, 1, 1, nz, ny, nx, nT) |>xpu;
ddvec= zeros(Float, nt, nr, 1, 1, 1, 1) |>xpu;

G!(d,m)= G_ps!(d, m, delay, nt, nr, nz, ny, nx, nT, wav0, mvec, dvec, Int(abs(lb)));
Gt!(grad, ∇d)= Gt_ps!(grad, ∇d, delay, nt, nr, nz, ny, nx, nT, wav0, g2, grvec, ddvec, Int(abs(lb)));

# eq_list=["hnd1", "cal1", "fij4", "nic", "fij5", "okt3", "chl1", "cps", "okt4", "rat", "spn", "chl2", "dnl", "kam", "fij1", "cal3", "nbl", 
#     "pb2", "crt", "van1", "chl3", "okt6", "cal2", "okt5", "mnd", "cal4", "nzd", "fij3", "csbg", "fij2", "okt1", "okt2", "bon1"];


heatmap(1:nr, tgrid, cpu(dobs), yflip= true, title= "dobs")

m2= zeros(length(mvec)) |>xpu;

λ= 0.001;
fill!(dg, 0.);
fill!(img, 0.);
algo= LBFGS(linesearch= BackTracking(order= 3))
backproject!(G!, Gt!, img, dobs, λ, mvec, m2, dg, algo, 5); 

fname= "results/"*eq*"_inv"
f= h5open(fname*".h5", "w")
write(f, eq, cpu(img));
close(f);