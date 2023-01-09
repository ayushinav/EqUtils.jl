using Plots
using HDF5
Float= Float64
include("utils_bp.jl");

dx= 10; dy= 10; 
limx= 500; limy= 500; # km
dT= 0.5; dt= 0.5;
nt= 801;
Tgrid = range(-50, 200, step= dT);
nT= length(Tgrid)

zgrid = Float.([0])
ygrid = Float.(range(-limy, stop=limy, step=dy))
xgrid = Float.(range(-limx, stop=limx, step=dx))
mgrid= [zgrid, ygrid, xgrid, Tgrid];
# rgrid= [rz, ry, rx];
nz, ny, nx, nT= length.(mgrid)
Tshift= Int.(round.(Tgrid./dt));

function get_gif(img)
    lim= maximum(extrema(img))
    anim= @animate for i in 1:nT
        heatmap(xgrid./111 .+cntr[2], ygrid./111 .+ cntr[1], img[1,:,:,i], yflip= false, title= string(Tgrid[i]), 
            clims= (0, lim), c=:seismic, aspect= 1, size= (200, 200),
            xlabel= "East ->", ylabel= "North ->"
        )
    end
    # images come out right with yflip=true, implies algorithm is working fine. We need yflip= false to get north pointing upwards (ygrid= -limy: dy: limy)
    return anim
end

function get_timesteps(img, tstamps)
    lim= maximum(extrema(img));
    p= [];
        for i in tstamps
            p1= heatmap(xgrid./111 .+cntr[2], ygrid./111 .+ cntr[1], img[1,:,:,i], yflip= false, title= string(Tgrid[i]), 
            clims= (-lim, lim), c=:seismic, cbar= false, yticks= false);
            # scatter!(p1, [cntr[2]],[cntr[1]], ms=0, grid=true, label= false)
            # xlabel= "East ->", ylabel= "North ->"); #, size= (400,400));
            if i== tstamps[1] || i== tstamps[length(tstamps)÷2+1]
                plot!(p1, yticks= true)
            end
            push!(p, p1);
        end
    
    pc= scatter([0,0], [0,1], xlims= (1,1.1), clims= (-lim, lim), 
        cbar= true, yshowaxis= false, xshowaxis= false, yticks= false, 
        xticks= false, label= false, c=:seismic, zcolor= [0,3]);
    push!(p, pc);
    peak_amp= [maximum(img[:,:,:,iT]) for iT in 1:nT];;
    pt= plot(Tgrid, peak_amp, label= "peak amplitude", xlabel= "time wrt origin (s)")
    push!(p, pt);
    layed= @layout [[grid(2,length(tstamps)÷2) a{0.1w}]
    b{0.3h}]
    plt= plot(p..., layout= layed, size= (1200, 800), guidefontsize= 10)
end

d0= 100000.
v0= 8.0;

for eqs in readdir("/home/abhinav/workspace/results/")[2:end]
    eq= split(eqs, ".")[1];
    
    f= h5open("results/"*eq*".h5", "r")
    # Synthetics
    img= read(f["chl1"]);
    close(f)
    cntr, _, _, _, _, _, _= obspy_data("chl1", "/home/abhinav/workspace/Bp_eq_data_processed", "/mnt/data2/isha/unique_pwin33_traces/");
    
    tstamps= range(start= 71, step= 10, length= 8)
    p= plot_result(img, tstamps)
    savefig("view_results/"*eq*".png")
end