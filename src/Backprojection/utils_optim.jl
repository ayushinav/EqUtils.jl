
"""`loss_L2!(m, d, dobs, nr)` </br>
Returns ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`"""
function loss_L2!(m, d, dobs, nr)
    G!(d, m);
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    J= Flux.mse(reshape(d, :), reshape(dobs, :));
    return J
end

"""`gradient_L2!(grad, m, d, dobs, N)` </br>
Returns the gradient of ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`, `N`= `nz`*`ny`*`nx`*`nT`"""
function gradient_L2!(grad, m, d, dobs, N)
    G!(d, m);
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    
    broadcast!(-, d, d, dobs);        
    Gt!(grad, d);
    rmul!(grad, 0.5*inv(N))
end

# ============================== damped LSS ======================================

"""`loss_L2!(m, d, dobs, nr)` </br>
Returns ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`"""
function loss_L2!(m, d, dobs, λ, nr)
    G!(d, m);
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    J= Flux.mse(reshape(d, :), reshape(dobs, :))+ λ*inv(length(m[:]))* Flux.norm(m)^2;
    return J
end

"""`gradient_L2!(grad, m, d, dobs, N)` </br>
Returns the gradient of ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`, `N`= `nz`*`ny`*`nx`*`nT`"""
function gradient_L2!(grad, m, d, dobs, λ, N)
    G!(d, m);
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    
    broadcast!(-, d, d, dobs);        
    Gt!(grad, d);
    rmul!(grad, 0.5*inv(length(dobs[:])))
    broadcast!(+, grad, grad, broadcast(*, 2λ*inv(length(m[:])), m))
end

# =========================== filtered x =========================================

"""`loss_L2!(m, d, dobs, nr)` </br>
Returns ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`"""
function loss_L2!(m, d, filt, dobs, λ, nr)
    G!(d, filt.(m));
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    J= Flux.mse(reshape(d, :), reshape(dobs, :))+ λ*inv(length(N))* Flux.norm(filt.(m))^2;
    return J
end

"""`gradient_L2!(grad, m, d, dobs, N)` </br>
Returns the gradient of ``||d- d_{obs}||^2``
where `d` will be updated inside the function as `d`= G`m`, `N`= `nz`*`ny`*`nx`*`nT`"""
function gradient_L2!(grad, m, d, filt, dfilt, dobs, λ, N)
    G!(d, filt.(m));
    # for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
    
    broadcast!(-, d, d, dobs);        
    Gt!(grad, d);
    rmul!(grad, 0.5*inv(N))
    broadcast!(+, grad, grad, broadcast(*, 2λ*inv(length(m[:])), m))
    broadcast!(*, grad, grad, dfilt.(m))
    # copyto!(grad, filt.)
end

# =========================== Frequency domain functions ==========================

# function func(x, ŷ, nt, nr)
#     Re= view(x, 1:nt*nr)
#     Im= view(x, nt*nr+1:2*nt*nr)
    
#     Re2= broadcast(*, Re, Re)
#     Im2= broadcast(*, Im, Im)
    
#     Norm2= broadcast(+, Re2, Im2, 0.001)
#     Norm= broadcast(sqrt, Norm2)
#     J= Flux.mse(Norm, ŷ)
#     return J;
# end

# A(x)= func(x, abs.(fobs[:]), nt, nr);
# dA(x)= Flux.gradient(A, x)[1];

# """`lossf_L2!(m, d, nr, H)` </br>
# Returns ``|| |H*d|- |H*d_{obs}| ||^2``
# where `d` will be updated inside the function as `d`= `H`G`m` and `H` is the Fourier matrix."""
# function lossf_L2!(m, d, nr, H)
#     G!(d, m);
#     for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
#     d_f= H*d;
#     x= cat(real.(dfc[:]), imag.(dfc[:]), dims= 1);
#     J= A(x);
#     return J
# end

# """`lossf_L2!(m, d, nr, H)` </br>
# Returns the gradient of ``|| |H*d|- |H*d_{obs}| ||^2``
# where `d` will be updated inside the function as `d`= `H`G`m` and `H` is the Fourier matrix."""
# function gradientf_L2!(grad, m, d, N)
#     G!(d, m);
#     for i in 1: nr rmul!(view(d,:,i), inv(0.001+std(view(d,:,i)))) end
#     dfc= H*dc;
#     x= cat(real.(dfc[:]), imag.(dfc[:]), dims= 1);
#     gr= dA(x)[1];
#     g= broadcast(+, view(x, 1:nt*nr), broadcast(*, im, view(x, nt*nr+1:2*nt*nr)));
#     g= reshape(g, nt,nr)
#     Gt!(grad, real.(H'*g));
# end
