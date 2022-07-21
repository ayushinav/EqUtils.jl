struct MyCustomException <: Exception 
    var::String
end

"""
Load the utitlies (center, θ, ϕ, PREM travel times) for `earthquake`

Right now, this has been customised for a particular file location. Since loading takes time, the utilities are then stored as HDF5 files at a local location for future use.

Returns center, θ, ϕ, PREM travel times, and network IDs for all the receivers. </br>
Note: θ, ϕ, and PREM travel times are with respect to the center
"""
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
                
        close(f)
        @info("Loading the grid data for "*earthquake*" from the loaded dataset in the saved directory.")
        return cntr, θ, ϕ, travel_tines, network
    
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
                f= open("earthquake.txt", "w") # tmp file to let python know which earthquake needs to be retrieved
                write(f, earthquake)
                close(f)

                run(`python3 read_obspy_data.py`);

                f= h5open(dir_path*"/"*earthquake*".h5")
                cntr= read(f, "center");
                θ= read(f, "theta")
                ϕ= read(f, "phi")
                travel_tines= read(f, "travel_times")
                network= read(f, "network")
                network= string.(network)
                close(f)
                @info(("Utilities for "*earthquake*" loaded and stored for future use."))
            end
        else
            throw(MyCustomException("The given earthquake "*earthquake*" does not exist in the data loadable directory."))
            return [];
        end
    end
    
    return cntr, θ, ϕ, travel_times, network;
end

"""
Load the traces for `earthquake`

Right now, this has been customised for a particular file location. Since loading takes time, the utilities are then stored as HDF5 files at a local location for future use.
"""

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

"""
Filters the data for particular arrays `arr`, and also changes the values of other utilities.

Arguments </br>
`arr`= List of stations, eg, ["TA", "AU"] </br>
`cntr`= Center of source grid, remains unchanged </br>
`θ`, `ϕ`, `travel_times`= utils that will be updated </br>
`network`= List of network IDs that will be updated </br>


Returns `cntr`, updated `θ`, updated `ϕ`, updated `travel_times`, updated list of `network`s, new number of receivers, `IDs` </br>
`IDs`= indices required to filter the same data from observed data later
"""
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


"""
Filters the data so as to get one receiver each from every array, and also changes the values of other utilities.

Arguments </br>
`cntr`= Center of source grid, remains unchanged </br>
`θ`, `ϕ`, `travel_times`= utils that will be updated </br>
`network`= List of network IDs that will be updated </br>

Returns `cntr`, updated `θ`, updated `ϕ`, updated `travel_times`, updated list of `network`s, new number of receivers, `IDs` </br>
`IDs`= indices required to filter the same data from observed data later
"""
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


"""
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
```
"""
function get_unique_receivers_per_bin(cntr, θ, ϕ, travel_times, network, Angles, bins, nr)
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
    @info("Number of receivers reduced to "*string(length(IDs)))
    return cntr, θ[IDs], ϕ[IDs], travel_times[IDs], network[IDs], length(IDs), IDs;
end

# =================================== Shift array ============================================


"""
Get shifts for all the receivers for all the points using Fraunhofer approximation.

Arguments </br>
`mgrid`= Bundle of source grid [`zgrid`, `ygrid`, `xgrid`, `Tgrid`] </br>
`rgrid`= Bundle of [`θ`, `ϕ`] </br>
`v0`= Velocity around the source grid **(Would be an important hyperparameter)**

Returns `shift` matrix of size (nz, ny, nx, nr)
"""
function get_shifts(mgrid, rgrid, v0)
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


# ========================== Cartesian- Polar conversion ========================================

Re= 6371; # radius of earth in km

"""
Converts from [z,y,x] to [θ, ϕ, depth] for `Re`= 6371 km.
"""
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

"""
Converts from [θ, ϕ, depth] to [z,y,x] for `Re`= 6371 km.
"""
function polar2cartesian(polar)
    θ, ϕ, dep= polar;
    r= Re- dep
    z= r*cos(π/180 *θ)
    y= r*sin(π/180 *θ)*sin(π/180 *ϕ)
    x= r*sin(π/180 *θ)*cos(π/180 *ϕ)
    return z,y,x
end