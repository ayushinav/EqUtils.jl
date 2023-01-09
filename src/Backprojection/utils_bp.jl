struct MyCustomException <: Exception 
    var::String
end

"""`obspy_data(earthquake, dir_path, load_path)` </br>
Load the utitlies (center, θ, ϕ, PREM travel times) for `earthquake`

Since loading takes time, the utilities are then stored as HDF5 files at a local location `dir_path` for future use. If the data is being loaded for the first time, the path of the earthquake data
`load_path` should be provided.

Returns center, θ, ϕ, PREM travel times, and network IDs for all the receivers. </br>
Note: θ, ϕ, and PREM travel times are with respect to the center"""
function obspy_data(earthquake, dir_path, load_path)
    eq_dir= readdir(dir_path);
    if earthquake*".h5" ∈ eq_dir
        f= h5open(dir_path*"/"*earthquake*".h5")
        cntr= read(f, "center");
        θ= read(f, "theta")
        ϕ= read(f, "phi")
        travel_tines= read(f, "travel_times")
        network= read(f, "network")
        network= string.(network)
        ids= read(f, "ids")
        ids= string.(ids);
        nevids=Int64.(read(f, "nevids"));
        close(f)
        @info("Loading the grid data for "*earthquake*" from the loaded dataset in the saved directory.")
        return cntr, θ, ϕ, travel_tines, network, ids, nevids
    
    else
        @info("Data for "*earthquake*" not pre-loaded.\n Loading it from loadable data directory.")

        # Not required now
        # eq_list= readdir(load_path) 
        # if false && earthquake ∈ eq_list
        #     eql= readdir(load_path*"/"*earthquake)
        #     trace_idx= -1
        #     for i in 1:length(eql)
        #         if eql[i][end-1:end]== ".a"
        #             trace_idx= i;
        #             break;
        #         end
        #     end
        #     if trace_idx== -1
        #         throw(MyCustomException(".a file for this earthquake does not exist."))
        #     else
                
        f= open("earthquake.txt", "w") # tmp file to let python know which earthquake needs to be retrieved
        write(f, earthquake)
        close(f)

        run(`python3 workspace/read_obspy_data.py`);

        f= h5open(dir_path*"/"*earthquake*".h5")
        cntr= read(f, "center");
        θ= read(f, "theta")
        ϕ= read(f, "phi")
        travel_tines= read(f, "travel_times")
        network= read(f, "network")
        network= string.(network)
        ids= read(f, "ids")
        ids= string.(ids);
        nevids= Int64.(read(f, "nevids"));
        close(f)
        @info(("Utilities for "*earthquake*" loaded and stored for future use."))
            # end
        # else
        #     throw(MyCustomException("The given earthquake "*earthquake*" does not exist in the data loadable directory."))
        #     return [];
        # end
    end
    
    return cntr, θ, ϕ, travel_times, network, ids. nevids;
end

"""`obspy_traces(earthquake, dir_path, load_path)` </br>
Load the traces for `earthquake`

Since loading takes time, the utilities are then stored as HDF5 files at a local location `dir_path` for future use. If the data is being loaded for the first time, the path of the earthquake data
`load_path` should be provided."""

function obspy_traces(earthquake, dir_path, load_path)
    eq_dir= readdir(dir_path);
    if earthquake*".h5" ∈ eq_dir
        f= h5open(dir_path*"/"*earthquake*".h5")
        d_obs= read(f, "traces");
        close(f)
        @info("Loading the data for "*earthquake*" from the loaded dataset in the saved directory.")
        return d_obs
    
    else
        @info("Data for "*earthquake*" not pre-loaded.\n Loading it from loadable data directory.")
        eq_list= readdir(load_path)
        if earthquake ∈ eq_list
            eql= readdir(load_path*"/"*earthquake)
            trace_idx= -1
            for i in 1:length(eql)
                if eql[i][end-1:end]== ".a"
                    trace_idx= i;
                    break;
                end
            end
            if trace_idx== -1
                throw(MyCustomException(".a file for this earthquake does not exist."))
            else
                f= open("earthquake.txt", "w")
                write(f, earthquake)
                close(f)

                run(`python3 read_obspy_traces.py`);

                f= h5open(dir_path*"/"*earthquake*".h5")
                d_obs= read(f, "traces");
                close(f)
                @info(("Processed traces for "*earthquake*" loaded and stored for future use."))
                return d_obs;
            end
        else
            throw(MyCustomException("The given earthquake "*earthquake*" does not exist in the loadable data directory."))
            return [];
        end
    end
end

# ================================== Filter arrays ===============================================

"""`filter_arrays(arr, cntr, θ, ϕ, travel_times, network, nr)` </br>
Filters the data for particular arrays `arr`, and also changes the values of other utilities.

Arguments </br>
`arr`= List of stations, eg, ["TA", "AU"] </br>
`cntr`= Center of source grid, remains unchanged </br>
`θ`, `ϕ`, `travel_times`= utils that will be updated </br>
`network`= List of network IDs that will be updated </br>


Returns `cntr`, updated `θ`, updated `ϕ`, updated `travel_times`, updated list of `network`s, new number of receivers, `IDs` </br>
`IDs`= indices required to filter the same data from observed data later"""
function filter_arrays(arr, cntr, θ, ϕ, travel_times, network, nr)
    IDs= [];
    for i in 1:length(arr)
        idx= (network.==arr[i]);
        ids= (1:nr)[idx];
        append!(IDs, ids)
    end
    @info("Number of receivers reduced to "*string(length(IDs)))
    return cntr, θ[IDs], ϕ[IDs], travel_times[IDs], network[IDs], length(IDs), IDs
end


"""`get_unique_receivers(cntr, θ, ϕ, travel_times, network, nr)` </br>
Filters the data so as to get one receiver each from every array, and also changes the values of other utilities.

Arguments </br>
`cntr`= Center of source grid, remains unchanged </br>
`θ`, `ϕ`, `travel_times`= utils that will be updated </br>
`network`= List of network IDs that will be updated </br>

Returns `cntr`, updated `θ`, updated `ϕ`, updated `travel_times`, updated list of `network`s, new number of receivers, `IDs` </br>
`IDs`= indices required to filter the same data from observed data later"""
function get_unique_receivers(cntr, θ, ϕ, travel_times, network, nr)
    IDs= [];
    for i in unique(network)
        iids= (1:length(network))[network.== i]
        iid= iids[1];
        # print(iid)
        append!(IDs, iid)
    end
    @info("Number of receivers reduced to "*string(length(IDs)))
    return cntr, θ[IDs], ϕ[IDs], travel_times[IDs], network[IDs], length(IDs), IDs;
end


"""`get_unique_receivers_per_bin(cntr, θ, ϕ, travel_times, network, Angles, bins, nr)` </br>
Filters the data so as to get one receiver each from every azimuthal bin, and also changes the values of other utilities.

Arguments </br>
`cntr`= Center of source grid, remains unchanged </br>
`θ`, `ϕ`, `travel_times`= utils that will be updated </br>
`network`= List of network IDs that will be updated </br>
`Angles`= Azimuths (not the same as of earthquakes) </br>
`bins`= Bins for which the uniform distribution is desired (One receiver per bin) </br>

Returns `cntr`, updated `θ`, updated `ϕ`, updated `travel_times`, updated list of `network`s, new number of receivers, `IDs` </br>
`IDs`= indices required to filter the same data from observed data later

`Angles` can be calculated through
```
γ= [[sin(θ[i]), sin(ϕ[i])*cos(θ[i]), cos(ϕ[i])*cos(θ[i])] for i in 1:nr]; # Take-off angles

Angles= zeros(nr)
for i in 1:nr
    γi= γ[i];
    angle= acos(γi⋅[1,0,0]);
    if γi[3]<0
        angle-= π;
    end
    Angles[i]= angle;
end
```"""
function get_unique_receivers_per_bin(cntr, θ, ϕ, travel_times, network, Angles, bins, nr, pert)
    IDs= [];
    bin_ids= zeros(nr)
    for ir in 1:nr
        bin_ids[ir]= bins[(Angles[ir].>=bins) .& (Angles[ir].< (bins.+step(bins)))][1];
    end
        
    for i in unique(bin_ids)
        iids= (1:length(bin_ids))[bin_ids.== i]
        iid= iids[1];
        # print(iid)
        append!(IDs, iid)
    end
    
    app= randn(length(IDs))
    app= Int.(round.(app./maximum(app).*pert));
    IDs= IDs+ app;
    IDs[IDs.<1].= 1;
    @info("Number of receivers reduced to "*string(length(IDs)))
    return cntr, θ[IDs], ϕ[IDs], travel_times[IDs], network[IDs], length(IDs), IDs;
end


function filter_nevids(eq, θ, ϕ, travel_times, network, ids, nr, nevids)
    bad_nevids= Dict([("fij1", [27]),
        ("hnd1", [17,28,29]),
        ("okt5", [7,14]),
        ("bon1", [47,70]),
        ("rat", [16,17,22]),
        ("okt1", [7,14,15,18,19,32]),
        ("nbl", [8,15])
        ]);
    
    # bad_nevids= Dict()
    # k_prev, _, bbin= split(bad[1], " ");
    # bad_nevids[k_prev]= [bbin];
    # for i in 2:length(bad)
    #     k, _, bbin= split(bad[i], " ");
    #     if k==k_prev push!(bad_nevids[k_prev], bbin);
    #     else bad_nevids[k]= [bbin]; k_prev= k;
    #     end
    #     # push!(a[k], bbin)
    # end
    
    IDs= [];
    if eq in keys(bad_nevids)
        append!(bad_nevids[eq], [1,2,3])
        append!(IDs, (1:nr)[[!(nevids[ir] in parse.(Int32, bad_nevids[eq])) for ir in 1:nr]]);
    else
        append!(IDs, (1:nr)[[!(nevids[ir] in [1,2,3]) for ir in 1:nr]]);
    end
    @info("Number of receivers reduced to "*string(length(IDs)))
    return θ[IDs], ϕ[IDs], travel_times[IDs], network[IDs], ids[IDs], length(IDs), nevids[IDs], IDs;
end
    
# =================================== Shift array ============================================


"""`get_shifts_fraunhofer(mgrid, rgrid, v0)` </br>
Get shifts for all the receivers for all the points using Fraunhofer approximation.

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**

Returns `shift` matrix of size (nz, ny, nx, nr)"""
function get_shifts_fraunhofer(mgrid, rgrid, v0)
    zgrid, ygrid, xgrid, _= mgrid;
    θ, ϕ= rgrid;
    nz, ny, nx= length.([zgrid, ygrid, xgrid]);
    nr= length(θ);
    
    γ= [[sin(θ[i]), sin(ϕ[i])*cos(θ[i]), cos(ϕ[i])*cos(θ[i])] for i in 1:nr]; # Take-off angles
    ξ= [[z,y,x]- [cz,cy,cx] for z in zgrid, y in ygrid, x in xgrid]; # Vector to the source point from the center
    ψ= [ξ[iz, iy, ix]⋅ γ[ir] for iz in 1:nz, iy in 1:ny, ix in 1:nx, ir in 1:nr];  # Fraunhofer path difference
    
    # Fraunhofer shifts for the source grid
    δt0= [ψ[iz,iy,ix,ir]/v0 for iz in 1:nz, iy in 1:ny, ix in 1:nx, ir in 1:nr]
    for i in 1:nr
        δt0[:,:,:,i]= δt0[:,:,:,i];
    end

    if !(minimum(δt0)+ Tgrid[1] > tgrid[1] && maximum(δt0)+ Tgrid[end] < tgrid[end]) # check if the Fraunhofer shifts along with the temporal shifts lie in the tgrid
        @warn("Fraunhofer shifts and temporal shifts exceed the tgrid")
    end

    shift= Int.(round.((δt0./dt)));
    return shift;
end



"""`get_shifts(mgrid, rgrid, v0)` </br>
Get shifts for all the receivers for all the points using Fraunhofer approximation.

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**

Returns `shift` matrix of size (nz, ny, nx, nr)"""
function get_shifts(mgrid, rgrid, v0, d0)
    zgrid, ygrid, xgrid, _= mgrid;
    rz, ry, rx= rgrid;
    nz, ny, nx= length.([zgrid, ygrid, xgrid]);
    nr= length(rx);
    
    ξ= [[z,y,x] for z in zgrid, y in ygrid, x in xgrid]; # Vector to the source point from the center
    ψ= [norm(ξ[iz, iy, ix]-  [rz[ir], ry[ir], rx[ir]]) for iz in 1:nz, iy in 1:ny, ix in 1:nx, ir in 1:nr];  # Eucledian distance
    t0= d0/v0;
    # shifts for the source grid
    δt0= [ψ[iz,iy,ix,ir]/v0 for iz in 1:nz, iy in 1:ny, ix in 1:nx, ir in 1:nr]
    for i in 1:nr
        δt0[:,:,:,i]= δt0[:,:,:,i];
    end

    δt0.= δt0.- t0;
    shift= Int.(round.((δt0./dt)));
    return shift;
end

function get_symae_data(eq)
    # f= h5open("/mnt/data2/isha/files/virtual_each.hdf5");
    f= h5open("/mnt/data2/isha/files/virtual_each_pb10s.hdf5");
    @show sizeof(f)   
    dat= read(f, eq);
    k = collect(keys(dat));
    sort!(k);
    nr= length(k);

    dobs= [dat[ik] for ik in k];
    dobs= reduce(hcat, dobs) |>xpu;
    @show sizeof(dat), sizeof(dobs)
    return dobs, k;
end

# ========================== Cartesian- Polar conversion ========================================

Re= 6371; # radius of earth in km

"""`cartesian2polar(z,y,x)` </br>
Converts from [z,y,x] to [θ, ϕ, depth] for `Re`= 6371 km."""
function cartesian2polar(z,y,x) 
    r= sqrt(x*x+y*y+z*z)
    lat= 180/π * asin(z/r)
    long= 180/π * acos(x/(r*cos(π/180 *lat)))
    if y== 0
        long= 90.
    else
        long= y/abs(y)* long;
    end
    dep= (Re - r)/1000
    return lat, long, dep
end

"""`polar2cartesian(polar)` </br>
Converts from [θ, ϕ, depth] to [z,y,x] for `Re`= 6371 km."""
function polar2cartesian(polar)
    θ, ϕ, dep= polar;
    r= Re- dep
    z= r*cos(π/180 *θ)
    y= r*sin(π/180 *θ)*sin(π/180 *ϕ)
    x= r*sin(π/180 *θ)*cos(π/180 *ϕ)
    return z,y,x
end

"""`(θ, ϕ)` to the coordinate system we use"""
function to_zyx(θ, ϕ)
    θ, ϕ= [θ, ϕ].*π/180;
    return [cos(ϕ), sin(ϕ)*(sin(θ)), sin(ϕ)*(cos(θ))];
end