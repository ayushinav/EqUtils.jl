# Parallel Stencils
@parallel_indices (it, ir, iT) function G_kernel!(dd, mm, ddelay, rrick0, lb)
    # if rrick0[it+lb-400- ddelay[ir, iT]]!=0
        dd[it, ir, iT] += rrick0[it+lb-400- ddelay[ir, iT]]* mm[iT];
    # end
    return
end

"""`G_ps!(d, m, delay, nt, nr, nz, ny, nx, nT, rick0, mvec, dvec, lb)` </br>
Forward model function. Updates `d` with G`m`. Other variables passed are for performance."""
function G_ps!(d, m, delay, nt, nr, nz, ny, nx, nT, rick0, mvec, dvec, lb)
    fill!(dvec, Float(0));
    copyto!(mvec, m);
    for iz in 1:nz, iy in 1:ny, ix in 1:nx
        delayv= view(delay, 1, :, iz, iy, ix, :);
        mv= view(mvec, 1, 1, iz, iy, ix, :);
        dv= view(dvec, :, :, 1, 1, 1, :);
        @parallel (1:nt, 1:nr, 1:nT) G_kernel!(dv, mv, delayv, rick0, lb);
    end
    copyto!(d, dropdims((sum(reshape(dvec, nt,nr,nT), dims=3)), dims=(3)));
end


@parallel_indices (it, ir, iT) function Gt_kernel!(grad, dd, ddelay, rrick0, lb)
    # if rrick0[it+lb-400- ddelay[ir, iT]]!=0 
        grad[it, ir, iT]= rrick0[it+lb-400- ddelay[ir, iT]]* dd[it,ir];
    # end
    return
end


"""`Gt_ps!(grad, ∇d, delay, nt, nr, nz, ny, nx, nT, rick0, g2, grvec, ddvec, lb)` </br>
Updates `grad` with ``G^T`` `m`. Other variables passed are for performance."""
function Gt_ps!(grad, ∇d, delay, nt, nr, nz, ny, nx, nT, rick0, g2, grvec, ddvec, lb)
    fill!(grvec, Float(0.));
    copyto!(ddvec, ∇d);
    dv= view(ddvec, :, :, 1, 1, 1, 1);
    for iz in 1:nz, iy in 1:ny, ix in 1:nx
        delayv= view(delay, 1, :, iz, iy, ix, :);
        gv= view(grvec, 1, 1, iz, iy, ix, :);
        @parallel (1:nt, 1:nr, 1:nT) Gt_kernel!(g2, dv, delayv, rick0, lb);
        g3= dropdims((sum(g2, dims=[1,2])), dims=(1,2));
        broadcast!(+, gv, gv, g3);
    end
    copyto!(grad, grvec);
end

function backproject!(G!, Gt!, img, dobs, mvec, dg, algo, itr)
    Gt!(img, dobs);
    if itr==0 return end
    fill!(mvec, 0.);

    J(m)= loss_L2!(m, dg, dobs, nr);
    g!(grad, m)= gradient_L2!(grad, m, dg, dobs, nz*ny*nx*nT);
    
    opt= optimize(J, g!, mvec[:], algo, Optim.Options(iterations= itr, show_trace= true))
    copyto!(img, opt.minimizer);
    return;
end

function backproject!(G!, Gt!, img, dobs, λ, mvec, dg, algo, itr)
    Gt!(img, dobs);
    if itr==0 return end
    fill!(mvec, 0.);

    J(m)= loss_L2!(m, dg, dobs, λ, nr);
    g!(grad, m)= gradient_L2!(grad, m, dg, dobs, λ, nz*ny*nx*nT);
    
    opt= optimize(J, g!, mvec[:], algo, Optim.Options(iterations= itr, show_trace= true))
    copyto!(img, opt.minimizer);
    return;
end

function backproject!(G!, Gt!, img, dobs, filt, dfilt, λ, mvec, dg, algo, itr)
    Gt!(img, dobs);
    if itr==0 return end
    # fill!(mvec, 0.);

    J(m)= loss_L2!(m, dg, filt, dobs, λ, nr);
    g!(grad, m)= gradient_L2!(grad, m, dg, filt, dfilt, dobs, λ, nz*ny*nx*nT);
    
    opt= optimize(J, g!, mvec[:], algo, Optim.Options(iterations= itr, show_trace= true))
    copyto!(img, opt.minimizer);
    return;
end