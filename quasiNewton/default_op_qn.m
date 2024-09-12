function options = default_op_qn
% The function generates the default options structure for the quasiNewton
% procedure.
%   If an initial Hessian is provided, substitute it in the field
%   invHess0 (the default is the identity matrix):
    options.invHess0 = [];
%   The externalHessian field is a boolean. It indicates whether the inverse Hessian 
%   is provided by the user as a function of x (true), or estimated
%   iteratively by quasi-Newton methods without any previous knwledge (false).
    options.externalHessian = false;
%   The field finvHess is a function handle which computes the inverse
%   Hessian.
%   If options.externalHessian==true, then the input is the x vector
%   values.
%   If options.externalHessian==false, then the input is of the following
%   form: @(invHess, dx, dgrad), where the input variables are defined in the program
%   HessianBFGS.m.
    options.finvHess = @HessianBFGS;
%   The field f_termination is a function handle which determines the
%   termination condition of the optimization process. The input variables
%   have to be the same as those defined in the default termination
%   function (see terminate_reldif.m for the definition of the input variables):
    options.f_termination = @(dif_f, dif_x, fmin, gradmin, xmin) terminate_reldif(dif_f, dif_x, fmin, gradmin, xmin, 1e-8, 1e-8);
    % The maximal number of iterations in the main quasi-Newton algorithm:
    options.maxNiter = 1e5;
%   The following fields limit the number of iterations in the different
%   phases of the line-search procedure. They are used to terminate the 
%   line-search process when something looks wrong.
%   The field alpha_factor imposes a restriction on the length of the bracket
%   interval obtained in the bracketing phase. It represents the maximal allowed
%   ratio of the final alpha in the bracketing phase to the initial alpha.
%   alpha_factor_1st represents the maximal allowed ratio in the first
%   iteration, which may be desired to be larger.
%   The field Lsection_factor represents the minimal allowed ratio of the
%   final interval obtained in the bracketing phase and the initial bracket
%   which was obtained in the bracketing phase.
%   The field Lsection_factor_1st represents the minimal allowed ratio in
%   the first iteration, which may be desired to be smaller.
    options.alpha_factor_1st = 1e10;
    options.Lsection_factor_1st = 1e-15;
    options.alpha_factor = 1e5;
    options.Lsection_factor = 1e-5;
%   The fields ro, sigma, tau1, tau2, tau3, are line-search parameters which are
%   defined in "Practical Methods of Optimization" by Fletcher, Sec. 2.6.
    options.ro = 0.01;
    options.sigma = 0.9;
    options.tau1 = 9;
    options.tau2 = 0.1;
    options.tau3 = 0.5;
%   The field minimal_f is defined as \bar{f} in the reference above.
    options.minimal_f = -Inf;
%   The field f_max_alpha is a function handle. It respresents a
%   restriction imposed on the allowed values of x in the line-search. The common
%   restriction is a limit on the maximal allowed magnitude of the x
%   values.
%   In order to limit the magnitude of x in the line-search to max_x, set:
%   options.f_max_alpha = @(x0, direction) alpha_max_x(x0, direction, max_x)
%   Any other restriction has to be of the form @(x0, direction), where x0
%   and direction are defined in alpha_max_x.m.
%   The default is no restriction.
    options.f_max_alpha = [];
%   The field Deltaf0 represents an initial estimation of the magnitude of the difference in
%   f between adjacent iterations, if available. Used for the determination of
%   the initial alpha value. Note that Deltaf0 is a positive number.
%   The default is no initial estimation for Deltaf0.
    options.Deltaf0 = [];
%   The plot field is a function handle which represents a plotting
%   function. It is used when it is desired to plot the convergence curve
%   during the process in each iteration. The arguments are the same as the
%   MATLAB function plot.
%   In order to plot in each iteration with the plot function, set:
%   options.plot = @plot
%   The default is no plotting.   
    options.plot = [];
%   The field alpha_estimation is a boolean. It indicates whethter the next
%   alpha is estimated as in the reference above (true), or chosen as the
%   resulting alpha in the search of the last iteration (false).
    options.alpha_estimation = true;
%   The field alpha0 is an initial value of alpha in the first iteration.
%   This is an alternative to the estimation by Deltaf0. If both alpha0 and
%   Deltaf0 are empty, the default initial alpha is 1.
    options.alpha0 = [];
end