function options = optionsOCqn(tolx, maxNiter)
% The function generates the options structure for the quasiNewton
% procedure for simple optimal control problems.
% tolx: The tolerance of the relative difference of the field
% maxNiter: The maximal number of iterations in the iterative search
% process
    options.invHess0 = [];
    options.externalHessian = false;
    options.finvHess = @HessianBFGS;
    options.f_termination = @(dif_f, dif_x, fmin, gradmin, xmin) norm(dif_x)/norm(xmin)<=tolx;
    options.maxNiter = maxNiter;
    options.alpha_factor_1st = 1e7;
    options.Lsection_factor_1st = 1e-15;
    options.alpha_factor = 1e5;
    options.Lsection_factor = 1e-10;
    options.ro = 0;
    options.sigma = 0.9;
    options.tau1 = 9;
    options.tau2 = 0.1;
    options.tau3 = 0.5;
    options.minimal_f = -Inf;
    options.f_max_alpha = [];
    options.Deltaf0 = [];
    % The plot function plots the maximization curve for the objective J:
    options.plot = @(x, y) plot(x, -y);
    options.alpha_estimation = true;
    options.alpha0 = [];
    % For an approximation of the Hessian by the Hessian of the penalty
    % term, use:
%     options.externalHessian = true;
%     options.finvHess = @(x) (here put the inverse.)
end