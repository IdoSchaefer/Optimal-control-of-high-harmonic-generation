function invHess_new = HessianBFGS(invHess, dx, dgrad)
% The function computes the updated inverse Hessian in the BFGS method.
% invHess: The previous inverse Hessian
% dx: The difference of the new solution vector from the old one (delta in Fletcher)
% dgrad: The difference of the new gradient vector from the old one (g in Fletcher)
    dx_dot_dgrad = dx.'*dgrad;
    invHess_dgrad = invHess*dgrad;
    dx_scaledT = dx.'/dx_dot_dgrad;
    invHess_dgrad_dxsT = invHess_dgrad*dx_scaledT;
    invHess_new = invHess + ((1 + (dgrad.'*invHess_dgrad)/dx_dot_dgrad)*(dx*dx_scaledT) - (invHess_dgrad_dxsT + invHess_dgrad_dxsT.'));    
end