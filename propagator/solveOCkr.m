function [U, field, mniter, matvecs] = solveOCkr(ihfun, K, V, x, ui, tdomain, Nts, Nt_ts, Nkr, tol, varargin)
% The program solves Schreodinger equation for a time dependent Hamiltonian,
% for an optimal control problem.
% It is based on Hillel Tal-Ezer's semi-global propagator.
% The function of matrix is computed by the Arnoldi approach.
% Input:
% K: the p diagonal element of the unperturbed Hamiltonian (kinetic
% energy).
% V: the time-independent potential. It's a function handle of the form:
% @(x). It is possible to insert the vector of the potential itself instead.
% Vtfun: A function handle of the form: @(u, x, t, more arguments). Returns
% the x vector of the time dependent, nonlinear perturbation: Vt(u, x, t).
% ihfun: a function. It computes the inhomogeneous vectors, from the
% potential energy and the wave function in the sampling points.
% ui: the initial state vector.
% x: the x grid
% tdomain: the time interval of the propagation: [ti tf]
% Nts: the number of time steps. 
% Nt_ts: the number of interior Chebyshev time points in the time step, used during
% the computational process.
% Nkr: the number of vectors in the Krylov space, for coputation of the
% function of operator.
% tol: the desired tolerance of the convergence.
% Output:
% U: contains the solution in the points: ti:((tf - ti)/Nts):tf, in different columns.
% mniter: the mean number of iteration for a time step; should be close to
% 1 for ideal efficiency of the algorithm.
% matvecs: the number of Hamiltonian operations.
%%%% Important remark: This is an old version of the propagator. A newer
%%%% version with some improvements is available online as supplementary material for the paper:
%%%% "Semi-global approach for propagation of the time-dependent
%%%% Schreodinger equation for time-dependent and nonlinear problems".
%%%% The adjustment of the newer version to OCT problems is straitforward,
%%%% but cumbersome. This is still in my "to-do" list (on 28.1.20).
    tinit = tdomain(1);
    tf = tdomain(2);
    % The length of the time interval of the whole propagation(can be negative):
    T = tf - tinit;
    % The length of the time step interval:
    Tts = T/Nts;
    Nx = length(ui);
    if length(V) == 1
    % If V is a function handle:    
        Vvec = V(x);
    else
    % If V is a vector:
        Vvec = V;
    end
    U = zeros(Nx, Nts*(Nt_ts - 1) + 1);
    field = zeros(1, Nts);
    U(:, 1) = ui;
    v0kr = zeros(Nkr, 1);
    v0Dbase = zeros(Nkr, 1);
    % The Chebyshev points for expansion in time, in the domain in which the Chebyshev expansion
    % is defined: [-1 1]
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The Chebyshev points for expansion in time, in the domain of the time
    % variable (changing of variables):
    t_ts = 0.5*(tcheb+1)*Tts;
    % The interior time points of the current time step for interpolation,
    % and the next time step for extrapolation:
    t_2ts = [t_ts.'; Tts + t_ts(2:Nt_ts).'].';
    Vkr = zeros(Nx, Nkr + 1);
    Hkr = zeros(Nkr + 1, Nkr);
    P = zeros(Nkr, Nkr);
    D = zeros(Nkr, Nkr);
    % Computing the matrix of the Taylor polynomials in time.
    % timeMts contains the points in the current time step, and timeMnext
    % contains the points in the next time step:
    [timeMts, timeMnext] = maketM(t_2ts(2:end), Nt_ts);
    % Computing the coefficients of the transformation from the Newton
    % interpolation polynomial terms, to a Taylor like form:
    Cr2t = r2Taylor4(t_ts, Tts);
    % The solution in the interior points in the time step.
    % Every column represents an interior time point in the time step:
    Unew = zeros(Nx, Nt_ts);
    % The v_vecs vectors are defined in a recursive way, and contain information about
    % the time dependence of the s_ext vectors:
    v_vecs = zeros(Nx, Nt_ts+1);
    v_vecs(:, 1) = ui;
    % The 0'th order approximation is the first guess, for the first time step.
    % Every column represents an interior time point in the time step:
    Uguess = guess_ts1(ui, Nt_ts);
    allniter = 0;
    Niter = 100;
    for tsi = 1:Nts
%         % The time, represented by the interior time points:
%         t = tinit + Tts*(tsi - 1) + t_ts;
        % The first guess for the iterative process, for the convergence of the u
        % values. Every column represents an interior time point in the time step:
        Ulast = Uguess;
        Unew(:, 1) = Ulast(:, 1);
        % Starting an iterative process, until convergence:
        niter = 0;
        reldif = tol + 1;
        while (reldif>tol && niter<Niter)
            % Calculation of the inhomogeneous s_ext vectors:
            [s_ext, Vts, field(tsi)] = ihfun(Ulast, x, tsi, Vvec, varargin{:});
            % Calculation of the coefficients of the form of Taylor
            % expansion, from the coefficients of the Newton
            % interpolation in the points t_ts.
            % The divided differences are computed by the function divdif.
            % For numerical stability, we have to transform the time points
            % in the time step, to points in an interval of length 4:
            Cnewton = divdif(t_ts*4/Tts, s_ext);
            % Calculating the Taylor like coefficients:
            Ctaylor = zeros(Nx, Nt_ts);
            Ctaylor(:, 1) = Cnewton(:, 1);
            for Newtoni = 2:Nt_ts
                Ctaylor(:, 1:Newtoni) = Ctaylor(:, 1:Newtoni) + Cnewton(:, Newtoni)*Cr2t(Newtoni, 1:Newtoni);
            end
            % Calculation of the v_vecs vectors:
            for polyi = 2:(Nt_ts+1)
                v_vecs(:, polyi) = (-1i*Hpsi(K, Vts ,v_vecs(:, polyi-1)) + Ctaylor(:, polyi-1))/(polyi - 1);
            end
            eigE = get_eigE;
            % Calculation of the wave function in all the time points
            % within the time step:
            Unew(:, 2:Nt_ts) = UfromLamb(timeMts, t_ts(2:Nt_ts));
            % To check the convergence of u:
            reldif = norm(Unew(:, Nt_ts) - Ulast(:, Nt_ts))/norm(Ulast(:, Nt_ts));
            Ulast = Unew;
            niter = niter + 1;
        end
        if niter == Niter
            display('The program has failed to achieve the desired tolerance.')
            % In such a case, change Nts, Nt_ts and/or Nkr.
        end
        allniter = allniter + niter;
        U(:, ((tsi - 1)*(Nt_ts - 1) + 2):(tsi*(Nt_ts - 1) + 1)) = Unew(:, 2:Nt_ts);
        % The new guess is an extrapolation from the points within the
        % previous time step:
        Uguess(:, 1) = Unew(:, Nt_ts);
        Uguess(:, 2:Nt_ts) = UfromLamb(timeMnext, t_2ts((Nt_ts + 1):(2*Nt_ts - 1)));
        v_vecs(:, 1) = U(:, tsi*(Nt_ts - 1) + 1);
        if ~isfinite(U(1, tsi*(Nt_ts - 1) + 1))
            % It means the solution diverges.
            % In such a case, change Nts, Nt_ts and/or Ncheb.
            display('Error.');
            return
        end
    end    
    mniter = allniter/Nts;
    matvecs = allniter*(Nt_ts + Nkr);
    
    %%% Nested functions: %%%
    
    function eigE = get_eigE
        createHkr(v_vecs(:, Nt_ts + 1));
%         %%% Debug:
%         if ~min(min(isfinite(Hkr)))
%             keyboard
%         end
%         %%% End of debug 
        [P, D] = eig(Hkr(1:Nkr, :));
        v0kr(1) = norm(v_vecs(:, Nt_ts + 1));
        v0Dbase = P\v0kr;
        eigE = diag(D);
    end

    function U = UfromLamb(timeM, t_vec)
        U = v_vecs(:, 1:(end-1))*timeM + Vkr(:, 1:Nkr)*P*spdiags(v0Dbase, 0, Nkr, Nkr)*f_fun(eigE, t_vec, Nt_ts, tol); 
    end

    function createHkr(v0)
        Vkr(:, 1) = v0/norm(v0);
        for vj = 1:Nkr
            Vkr(:, vj+1) = ifft(K.*fft(Vkr(:, vj))) + Vts.*Vkr(:, vj);
            for vi = 1:vj
                Hkr(vi, vj) = Vkr(:, vi)'*Vkr(:, vj+1);
                Vkr(:, vj+1) = Vkr(:, vj+1) - Hkr(vi, vj)*Vkr(:, vi);
            end
            Hkr(vj+1, vj) = norm(Vkr(:, vj+1));
            Vkr(:, vj+1) = Vkr(:, vj+1)/Hkr(vj+1, vj);
        end
    end

end

%%% Sub functions: %%%

function result = f_fun(z, t, Nt_ts, tol)
% The function f(-1i*z, t):
    persistent factorialNt_ts
    if isempty(factorialNt_ts)
        % This is in order to avoid the repeated computation by the
        % factorial procedure, which is time-consuming:
        factorialNt_ts = factorial(Nt_ts);
    end
    Nt = length(t);
    Nz = length(z);
    minus_izt = -1i*z*t;
    % Condition for estimating if f_fun(z, t) should be computed directly or by
    % a "tail" of a Taylor expansion:
    is_big = factorialNt_ts*eps./abs(minus_izt.^(Nt_ts)) < tol;
    result = ones(Nz, Nt);
    % First, we compute f(-1i*z, t)/(t^Nt_ts), which is a function of zt.
    % A direct computation for large arguments:
    result(is_big) = exp(minus_izt(is_big));
    for polyi = 1:Nt_ts
        result(is_big) = polyi*(result(is_big) - 1)./minus_izt(is_big);
    end
    % Computation by a Taylor form for small arguments:
    is_not_converged = ~is_big;
    term = double(is_not_converged);
    polydeg = 1;
    while max(max(is_not_converged))
        term(is_not_converged) = minus_izt(is_not_converged).*term(is_not_converged)/(polydeg + Nt_ts);
        result(is_not_converged) = result(is_not_converged) + term(is_not_converged);
        polydeg = polydeg + 1;
        is_not_converged(is_not_converged) = abs(term(is_not_converged))./abs(result(is_not_converged)) > eps;
    end
    % Obtaining the required function f(-1i*z, t):
    result = result.*((ones(Nz, 1)*t).^Nt_ts);
end

function [timeMts, timeMnext] = maketM(t_2ts, Nt_ts)
% Computation of the matrix timeM, of time Taylor polynomials.
% timeM(i, j) = t_2ts(j)^(i - 1)
    Nt_2ts = 2*Nt_ts - 2;
    timeM = zeros(Nt_ts, Nt_2ts);
    timeM(1, :) = ones(1, Nt_2ts);
    for vi = 2:Nt_ts
        timeM(vi, :) = t_2ts.*timeM(vi - 1, :);
    end
    timeMts = timeM(:, 1:(Nt_ts - 1));
    timeMnext = timeM(:, Nt_ts:Nt_2ts);    
end

function Uguess = guess_ts1(ui, Nt_ts)
    Uguess = ui*ones(1, Nt_ts);
end