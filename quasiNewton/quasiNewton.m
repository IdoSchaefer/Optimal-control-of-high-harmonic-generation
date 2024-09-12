function [xmin, fmin, gradmin, niter, nfevals, dif_f, dif_x, conv, alpha_last, invHess] = quasiNewton(fDf, x0, options)
% The function minimizes a function which depends on a vector x by the BFGS method.
% The line-search is based on the procedure described in "Practical
% Methods of Optimization" by Fletcher, Sec. 2.6, for available gradient
% data.
% Input:
% fDf: A function handle which returns the minimized function and its
% gradient as two different output arguments in the following way:
% [function, gradient] = fDf(x), where the input x is a column vector.
% x0: The guess solution, a column vector 
% options: A structure of the optimization options. The program
% default_op_qn.m sets the defaults. An explanation of the various fields in the
% structure can be found therein.
% Output:
% xmin: The minimized solution
% fmin: The value of f at xmin
% gradmin: The gradient value at xmin
% niter: The number of required iterations
% nfevals: The number of f and fgrad evaluations
% dif_f: The difference of fmin from the f value at the previous iteration
% dif_x: The difference of xmin from the x value at the previous iteration
% (a vector).
% conv: A row vector of length niter + 1 which represents the convergence
% history; contains the f values in all iterations. The first entry is the
% value of f of the initial guess solution.
% alpha_last: The resulting alpha in the last iteration
% invHess: The resulting inverse Hessian in the last iteration
    % A boolean which indicates if the convergence curve is plotted in each
    % iteration:
    plot_mode = ~isempty(options.plot);
    niter = 0;
    % Initialization of the search from the guess solution:
    xmin_old =  x0;
    [fmin_old, gradmin_old] = fDf(x0);
    if nargout>7 || plot_mode
        % If the convergence information is required:
        conv = zeros(1, options.maxNiter + 1);
        conv(1) = fmin_old;
    end
    % The dimension of the problem:
    dim = length(x0);
    % Setting the initial direction of search:
    if isempty(options.invHess0)
        % If the initial inverse Hessian is not supplied, it is taken as
        % the identity matrix:
        invHess = eye(dim);
        direction = -gradmin_old;
    else
        % If the initial inverse Hessian is supplied:
        invHess = options.invHess0;
        direction = -invHess*gradmin_old;
    end
    % This initialization is required for the check of the convergence
    % condition in the first iteration:
    dif_f = Inf;
    dif_x = Inf*ones(dim, 1);
    nfevals = 1;
    % A boolean which indicates if the solution f is has achieved a minimal
    % value, which is defined by the user as a sufficient condition to be
    % considered as an acceptable solution:
    sol_achieved = false;
    % A boolean which indicated if the user has chosen to stop the process,
    % when the plot mode is applied:
    stop = false;
    if plot_mode
        figure
        % In the plot mode, it is possible to terminate the process interactively,
        % or to pause the process and enter a debugging mode:
        uicontrol('style', 'push', 'string', 'Stop', 'callback', @buttonStop, 'position', [0, 0, 60, 20]);
        uicontrol('style', 'push', 'string', 'Pause', 'callback', @EnterPrompt, 'position', [70, 0, 60, 20]);
    end
    % Initialization of alpha in various user options:
    if ~isempty(options.alpha0)
        % If alpha is supplied directly:
        alpha1 = options.alpha0;
    elseif ~isempty(options.Deltaf0)
        % If the estimated magnitude of the difference in f between adjacent
        % iterations is supplied, then alpha is estimated by a quadratic
        % approximation (see Fletcher), with an upper bound 1:
        alpha1 = min([1, -2*options.Deltaf0/(gradmin_old.'*direction)]);
    else
        % The default, when no external estimation is provided:
        alpha1 = 1;
    end
    if isempty(options.f_max_alpha)
        % If no upper bound for alpha is supplied:
        max_alpha = Inf;
    else
        % If an upper-bound for the allowed alpha in the line-search is
        % defined by the user, as a function of the solution x and the direction of search:
        max_alpha = options.f_max_alpha(x0, direction);
        if alpha1>max_alpha
            alpha1 = max_alpha;
            fprintf('\nWarning: The alpha1 value (%d) was limited by the user preferences (iteration No. 1).\n', alpha1)
        end
    end
    % Setting bounds for the two phases in the line-search procedure (see
    % default_op_qn.m):
    alpha_factor = options.alpha_factor_1st;
    Lsection_factor = options.Lsection_factor_1st;
    % Starting the iterative search proceudre:
    while ~options.f_termination(dif_f, dif_x, fmin_old, gradmin_old, xmin_old) && niter<options.maxNiter && ~sol_achieved &&...
            alpha1>0 && ~isnan(fmin_old) && ~stop
    % The iterative procedure continues while the ending conditions defined
    % by the user have not been achieved yet (conditions 1-3), and the
    % procedure can proceed (conditions 4, 5), and the process has not been
    % stopped by the user (condition 6).
%       Use this for debugging if something looks wrong:
%       if niter==8
%           keyboard
%       end
        %tic
        % The line-search procedure:
        [xmin, fmin, gradmin, sol_achieved, nfevals_ls, alpha_last] = line_search1(fDf, xmin_old, fmin_old, gradmin_old, direction,...
            alpha1, options.ro, options.sigma, options.tau1, options.tau2, options.tau3, options.minimal_f, max_alpha,...
            alpha_factor, Lsection_factor);
        nfevals = nfevals + nfevals_ls;
        dif_f = fmin - fmin_old;
        dif_x = xmin - xmin_old;
        if options.externalHessian
            % If an external procedure for the computation of the inverse Hessian is supplied: 
            invHess = options.finvHess(xmin);
        else
            % A quasi-Newton estimation for the inverse-Hessian:
            dgrad = gradmin - gradmin_old;
            invHess = options.finvHess(invHess, dif_x, dgrad);
        end
        % The new direction of search:
        direction = -invHess*gradmin;
        if options.alpha_estimation
            % If alpha is chosen by the user to be estimated in each
            % iteration by a quadratic interpolation:
            Deltaf = max([-dif_f, 10*eps]);
            alpha1 = min([1, -2*Deltaf/(gradmin.'*direction)]);
        else
            % The alpha value is chosen as the last value, with an upper
            % bound 1:
            alpha1 = min([1, alpha_last]);
        end
        if ~isempty(options.f_max_alpha)
            % If an upper-bound for the allowed alpha in the line-search is
            % defined by the user:
            max_alpha = options.f_max_alpha(xmin, direction);
            if alpha1>max_alpha
                alpha1 = max_alpha;
                fprintf('\nWarning: The alpha1 value (%d) was limited by the user preferences (iteration No. %d).\n', alpha1, niter + 2)
            end
        end
        xmin_old = xmin;
        fmin_old = fmin;
        gradmin_old = gradmin;
        niter = niter + 1;
        if niter==1
            alpha_factor = options.alpha_factor;
            Lsection_factor = options.Lsection_factor;
            % The values in the first iteration are typically chosen to be
            % different from the other iterations, and thus have to be
            % updated after the first iteration.
        end
        if nargout>7 || plot_mode
            % Updating the convergence information, if necessary:
            conv(niter + 1) = fmin;
        end
        if plot_mode
            % The updated convergence curve is plotted:
            options.plot(0:niter, real(conv(1:(niter + 1))))
            drawnow
        end
        %toc
    end
    if sol_achieved
        fprintf('\nThe optimization process has achieved the minimal function value.\n')
    elseif niter==options.maxNiter
        fprintf('\nOptimization failed. The maximal number of iterations has been achieved.\n')
    elseif max_alpha<=0
        fprintf('\nThe optimization process was stopped since the limitation on the magnitude of the alpha parameter does not allow further progress.\n')
    elseif alpha1<=0
        fprintf('\nThe optimization process was stopped since the direction of search is not a descent direction.\n')
    elseif isnan(fmin_old)
        fprintf('\nWarning: The function value is NaN or Inf. The process was stopped.\n')
    elseif stop
        fprintf('\nThe optimization process was stopped by the user.\n')
    end
    if nargout>7
        conv = conv(1:(niter + 1));
    end
    
    %%%% Nested function: %%%%

    function buttonStop(hObject, event)
        % Applied when the user presses the "Stop" button on the figure:
        stop = true;
    end
    
end

%%%% Sub function: %%%%

function EnterPrompt(hObject, event)
% Applied when the user presses the "Pause" button on the figure:
    keyboard
end