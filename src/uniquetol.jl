

function uniquetol_list(itr::Array{T,1},tol=1e-6) where {T<:Complex}
    preal = sortperm(real(itr))     # buffer for permutation idxs of reals
    dreal = diff(real(itr[preal]))  # buffer for differences of reals

    pimag = Array{Int}(size(preal))     # buffer for permutation idxs of imags
    dimag = Array{Float64}(size(dreal)) # buffer for differences of imags

    out = [ ];
    pri1=1;
    for pri2 = 1:length(preal)
        # look for tol-wide gap in reals
        if pri2==length(preal) || dreal[pri2] > tol
            pridxs = preal[pri1:pri2]
            sortperm!( pimag[1:length(pridxs)], imag(itr[pridxs]) )
            dimag[1:(length(pridxs)-1)] = diff(imag(itr[preal[pridxs]]))

            pii1=1;
            for pii2 = 1:length(pridxs)
                if pii2==length(pridxs) || dimag[pii2] > tol
                    push!(out, Array{T,1}(itr[pridxs[pii1:pii2]]))
                    pii1=pii2+1;
                end # imag tol check
            end # pii2, tol-large gap in imags

            pri1=pri2+1;
        end # real tol check
    end # pri2, tol-large gap in reals
    return out
end


function uniquetol_list(itr::Array{T,1},tol=1e-6) where {T<:Real}
    preal = sortperm(real(itr))     # buffer for permutation idxs of reals
    dreal = diff(real(itr[preal]))  # buffer for differences of reals

    out = [ ];
    pri1=1;
    for pri2 = 1:length(preal)
        # look for tol-large gap in reals
        if pri2==length(preal) || dreal[pri2] > tol
            push!(out, Array{T,1}(itr[preal[pri1:pri2]]))
            pri1=pri2+1;
        end # real tol check
    end # pri2, tol-large gap in reals
    return out
end

