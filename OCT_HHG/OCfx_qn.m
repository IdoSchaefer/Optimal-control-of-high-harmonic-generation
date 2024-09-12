function [fieldt, fieldw, psi, evat, evaw, evmut, evmuw, relE, conv, niter, mallniterc, Jmax, maxgrad, alpha, invHess] = OCfx_qn(psi0,...
    Vf, m, xdomain, mu_x, a0, penalnormf, Dpenalnormf, fguess, filterE, filtera, options, T, dt, Nt_ts, Nkr, tol, maxNiter)
% The function solves an OCT procedure for the optimization of HHG in one dimension, or
% other HG problems defined on a one-dimensional spatial grid. A
% quasi-Newton optimization procedure is applied.
%%% Input: %%%
% psi0: The initial state vector
% Vf: The stationary potential; a functon handle of the form @(x), where
% the input is the x grid. Alternatively, Vf can be a column vector of the
% potential values in x.
% m: The mass (in HHG, the electron mass, which is 1 in atomic units)
% xdomain: The boundaries of the x domain, where xdomain = [min_x, max_x]
% mu_x: The dipole function mu evaluated at the x grid points; a column 
% vector of the dimension of the x grid; mu_x is used in the Hamiltonian for the
% interaction term with the external field, -field(t)*mu_x. In HHG, insert 
% here the numerical x vector which its derivative decays to 0 near the 
% absorbing boundaries.
% a0: the stationary acceleration, evaluated at the x grid points; a
% column vector of the dimension of the x grid 
% penalnormf: The penalty function on the ionization; a function handle of
% the form @(sq_norm_psi), where the input is <psi(T)|psi(T)>
% Dpenalnormf: The derivative of the function represented by penalnormf; a
% function handle of the form @(sq_norm_psi), where the input is <psi(T)|psi(T)>
% fguess: The guess field function in the frequency domain; a function
% handle of the form @(w), where w is a column vector which represents the
% omega grid. Aternatively, fguess can be a column vector of the guess
% function evaluated at the omega grid points. If the guess does not
% satisfy the 0 boundary conditions of the temporal profile, the program
% constructs a constrained field with the 0 boundary conditions from the
% unconstrained field.
% filterE: The scaled filter function of the driving electric field spectrum; a
% function handle of the form @(w), where w is a column vector which
% represents the omega grid
% filtera: The filter function of the emission spectrum, represented by the
% stationary acceleration expectation spectrum; a function handle of the
% form @(w), where w is a column vector which represents the omega grid
% options: The options structure of the BFGS search procedure (see
% quasiNewton.m and default_op_qn.m); if is empty (substituted with []), the default is set by
% the procedure optionsOCqn.m, using the tol and maxNiter input parameters.
% T: The final time
% dt: The time-step
% Nt_ts: The number of Chebyshev points in the internal grid in each
% time-step, using the semi-global propagator.
% Nkr: The dimension of the Krylov approximation in the semi-global propagator
% tol: The tolerance of the relative difference of the driving field; 
% active only if options is empty, but has to be defined anyway (this
% defect is corrected in OCfx_qn1.m).
% maxNiter: The maximal number of iterations in the BFGS optimization
% procedure; active only if options is empty (this defect is corrected in OCfx_qn1.m).
%%% Output: %%%
% fieldt: The optimized time-dependent electric field, evaluated at the
% points of the main time-grid, t = 0:dt:T (not at the intermediate 
% Chebyshev point structure of the semi-global propagation); a row vector
% fieldw: The optimized electric field in the omega domain; a row vector
% psi: The state in all time-points of the process under the influence of the optimized
% field; evaluated at the main time-grid, where the different time-points
% are represented by separate columns.
% evat: The time-dependent expectation value of the stationary
% acceleration, evaluated at the main time-grid; a row vector
% evaw: evat in the omega domain; a row vector
% evmut: The time-dependent expectation value of the dipole, defined by
% mu_x; a row vector
% evmuw: evmut in the omega domain; a row vector
% relE: The relative difference of the field from the previous iteration at the
% end of the optimization process
% conv: The convergence history---J as a function of the iteration number;
% a row vector of dimension niter+1, where the 1'st entry represents the
% guess field
% niter: The total number of iterations
% mallniterc: The mean number of semi-global iterations in the whole
% optimization process; should be close to 1 for efficiency.
% Jmax: The Jmax value of the optimized field
% maxgrad: The infinity norm of the gradient at the end of the optimization
% process
% alpha: The alpha parameter at the end of the optimization process (see
% quasiNewton.m)
% invHess: The inverse Hessian at the end of the optimization process

    % The number of time-points:
    Nt = T/dt;
    % The distance between adjacent points in the omega grid:
    dw = pi/T;
    % Integration weights of the boundary including omega grid:
    integw = [dw/2; ones(Nt - 1, 1)*dw; dw/2];
    % Additional factor for conversions between the time and frequency grid
    % by the DCT transformation:
    dctfactor = T/(sqrt(Nt*pi));
    % The dimension of the problem, whcih is the number of x grid points:
    Nx = length(psi0);
    % The index of the last point in the semi-global time-grid, where the
    % internal Chebyshev points are also counted:
    allt_lasti = Nt*(Nt_ts - 1) + 1;
    % Constructing the x grid:
    min_x = xdomain(1);
    max_x = xdomain(2);
    xdlength = max_x - min_x;
    dx = xdlength/Nx;
    x = (min_x:dx:(max_x - dx)).';
    % Constructing the p (momentum) grid:
    p = (0:(2*pi/xdlength):(2*pi*(1/dx - 1/xdlength))).';
    p((Nx/2 + 1):Nx) = p((Nx/2 + 1):Nx) - 2*pi/dx;
    % The kinetic energy in the p domain:
    K = p.^2/(2*m);
    % Defining the complex conjugate potential for the propagation of chi:
    if length(Vf) == 1
        % If Vf is a function handle:
        conjVf = @(x) conj(Vf(x));
    else
        % If Vf is a vector:
        conjVf = conj(Vf);
    end
    % The stationary acceleration diagonal matrix:
    a0M = spdiags(a0, 0, Nx, Nx);
    % The tolerance parameter of the propagation procedure:
    tolprop = 1e-3*tol;
    % The sum of the mean number of semi-global iterations of all 
    % propagations performed during the optimization process (required for
    % the computation of mallniterc):
    summniterc = 0;
    % The time-dependent field in the semi-global propagation grid (including the internal
    % Chebyshev points in each time-step):
    allfield = zeros(allt_lasti, 1);
    % The time-dependent stationary acceleration expectation value in the
    % semi-global propagation grid:
    allevat = zeros(allt_lasti, 1);
    % allevat filtered by the filtera function:
    evafil = zeros(allt_lasti, 1);
    % allevat in the omega domain:
    evaw = zeros(1, Nt + 1);
    % The integrand of the Jmax term in the functional:
    Jmax_fun = zeros(Nt + 1, 1);
    % The term <chi(t)|mu|psi(t)> evaluated at the semi-global propagation
    % grid:
    chimiupsi = zeros(allt_lasti, 1);
    % The reversed Chebyshev points in each time-step in the [-1, 1] domain:
    tcheb = -cos(((1:Nt_ts) - 1)*pi/(Nt_ts-1));
    % The internal Chebyshev points of each time-step in the [0, dt] domain:
    t_ts = 0.5*(tcheb+1)*dt;
    % The omega grid:
    w = (0:pi/T:pi/dt).';
    % The vector of the driving field filter-function, evaluated at w:
    vfilterE = filterE(w);
    % The indices of the w which participate in the optimization, where
    % the driving field spectrum does not attain negligible values (nz is
    % non-zero):
    iEnz = find(vfilterE>=eps*max(vfilterE));
%    [iEnz, ~, vfilterEnz] = find(vfilterE);
%    iEnz = find(vfilterE>=tolprop*max(vfilterE)/10);
    % vfilterE in the points which participate in the optimization:
    vfilterEnz = vfilterE(iEnz);
    % The integration weights of the internal grid Chebyshev points in the
    % semi-global prapagation time grid:
    igweights = chebweights(Nt_ts, 1);
    igweights = [2*igweights(1), igweights(2:(Nt_ts - 1))];
    % The integration weights of the omega grid points which participate in
    % the optimization:
    integwnz = integw(iEnz);
    % The inverse cosine transform of the scaled filter function, evaluated
    % at t=0:
    dctfilterE0 = sqrt(2/pi)*sum(vfilterEnz.*integwnz);
    % Constructing a vector of cos(w*T) evaluated at all w values:
    coswT = ones(Nt + 1, 1);
    coswT(2:2:(Nt + 1)) = -1;
    % The inverse cosine transform of the scaled filter function, evaluated
    % at t=T:
    dctfilterET = sqrt(2/pi)*sum(vfilterEnz.*coswT(iEnz).*integwnz);
    % The determinant of the M matrix in the paper:
    deTdctfilterE  = dctfilterE0^2 - dctfilterET^2;
    % The driving field spectrum in all w:
    fieldw = zeros(Nt + 1, 1);
    if length(fguess) == 1
        % if fguess is a function handle:
        for wnzi = iEnz
            fieldw(wnzi) = fguess(w(wnzi));
        end
    else
        % if fguess is a vector:
        fieldw(iEnz) = fguess(iEnz);
    end
    % The driving field spectrum the omega values which participate in the
    % optimization:
    fieldwnz = fieldw(iEnz);
    % If the guess field doesn't satisfy the temporal boundary conditions, a new
    % guess field is constructed which satisfies the boundary conditions:
    if sqrt(2/pi)*abs(sum(fieldwnz.*integwnz)) > tolprop*1e-2 || sqrt(2/pi)*abs(sum(fieldwnz.*coswT(iEnz).*integwnz)) > tolprop*1e-2
        fieldwnz = getfieldwcon(fieldwnz);
    end
    % The filter function of the stationary acceleration expectation,
    % evaluated at the w grid points:
    vfiltera = filtera(w);
    % psi at the semi-global propagation grid:
    allpsi = zeros(Nx, allt_lasti);
    % a*psi at the semi-global propagation grid, where a is the stationay
    % acceleration operator:
    allapsi = zeros(Nx, allt_lasti);
    % The inhomogeneous term in the inhomogeneous Schreodinger equation for
    % the chi propagation:
    chiihterm = zeros(Nx, allt_lasti);
    % nprop counts the number of propagations:
    nprop = 0;
    if isempty(options)
        % Setting the default options:
        options = optionsOCqn(tol, maxNiter);
    end
    if isempty(options.invHess0)
        % Initializing the approximated Hessian to the default: 
        options.invHess0 = diag(vfilterEnz./(2*integwnz));
    end
    % Performing an optimization process of fieldwnz:
    [fieldwnz, ~, minusgrad, niter, ~, ~, dif_fieldw, minus_conv, alpha, invHess] = quasiNewton(@Jeval, fieldwnz, options);
    fieldw = fieldw.';
    maxgrad = max(abs(minusgrad));
    conv = -minus_conv;
    relE = norm(dif_fieldw)/norm(fieldwnz);
    niter
    nprop
    psi = allpsi(:, 1:(Nt_ts - 1):allt_lasti);
    fieldt = allfield(1:(Nt_ts - 1):allt_lasti).';
    evat = allevat(1:(Nt_ts - 1):allt_lasti).';
    allevmiut = real(evmiu(allpsi, mu_x));
    evmut = allevmiut(1:(Nt_ts - 1):allt_lasti);
    evmuw = dctIfrom_ig(allevmiut, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
    mallniterc = summniterc/nprop;
    Jmax = sum(Jmax_fun);
    beep
    
    %%% Nested functions: %%%
    
    function [minusJ, minusgrad] = Jeval(fieldwnz)
    % The optimization function
    % Input:
    % fieldwnz: The driving field spectrum at the omega points which
    % participate in the optimization
    % Output:
    % minusJ: The objective of the optimization, which is -J
    % minusgrad: The gradient of the optimization, which is minus the
    % gradient of J w.r.t. fieldwnz
        fieldw(iEnz) = fieldwnz;
        % Obtaining the the temporal driving field from its spectrum,
        % evaluated at the semi-global propagation grid:
        allfield = dctIintgrid(fieldw, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        % The forward propagation of psi:
        [allpsi, ~, mniterc] = solveOCkr(@ihfieldmiux, K, Vf, x, psi0, [0 T], Nt, Nt_ts, Nkr, tolprop, allfield.', mu_x);
        summniterc = summniterc + mniterc;
        % Computing the stationary acceleration expectation at the semi-global
        % propagation grid, for the chi propagation:
        allapsi = a0M*allpsi;
        % The time-dependent part of the acceleration operator is omitted,
        % since it contributes low frequencies only. This just increases the
        % background of the signal, and does not contain any interesting
        % information.
        for allti = 1:allt_lasti
        % The values of the expectation value are supposed to be real:
            allevat(allti) = real(allpsi(:, allti)'*allapsi(:, allti));
        end
        % The spectrum of the stationary acceleration expectation is
        % computed utilizing the data in the whole propagation grid:
        evaw = dctIfrom_ig(allevat, T, t_ts(1:(Nt_ts - 1)), igweights)*dctfactor;
        % Obtaining chiihterm (the inhomogeneous term for the chi
        % propagation):
        get_chiihterm(allapsi);
        psiT = allpsi(:, allt_lasti);
        % The square norm of psi at t=T, for the computation of J:
        sq_norm_psiT = psiT'*psiT;
        % The final condition for the chi backward propagation:
        chiT = Dpenalnormf(sq_norm_psiT)*psiT;
        % The backward propagation of chi:
        [allchi, ~, mniterc] = solveOCkr(@ihalltchimiux, K, conjVf, x, chiT, [T 0], Nt, Nt_ts, Nkr, tolprop,...
            allfield(allt_lasti:-1:1).', mu_x, chiihterm);
        summniterc = summniterc + mniterc;
        nprop = nprop + 2;
        % Computing the objective -J:
        Jmax_fun = 0.5*evaw.^2.*vfiltera.*integw;
        Jenergy_fun = -(fieldwnz.^2./vfilterEnz).*integwnz;
        minusJ = -(sum(Jmax_fun) + sum(Jenergy_fun) + penalnormf(sq_norm_psiT));
        % Computing the -gradient.
        %  Obtaining <chi(t)|mu|psi(t)> at the propagation grid:
        for allti = 1:allt_lasti
            chimiupsi(allti) = -imag(allchi(:, allt_lasti - allti + 1)'*(mu_x.*allpsi(:, allti)));
        end
        % Computing chimiupsi transformed to the frequency domain,
        % utilizing the data in the entire propagation grid:
        dct_chimiupsi = dctfactor*dctIfrom_ig(chimiupsi, T, t_ts(1:(Nt_ts - 1)), igweights);
        % The unconstrained fieldw:
        newfieldw_unc_nz = dct_chimiupsi(iEnz).*vfilterEnz;
        % Computing the Lagrange-multiplies of the temporal boundary
        % constraints from the unconstrained field:
        [lambda0, lambdaT] = fieldw_unc2lambda(newfieldw_unc_nz);
        minusgrad = 2*integwnz.*(fieldwnz./vfilterEnz - dct_chimiupsi(iEnz) + lambda0 + lambdaT*coswT(iEnz));
    end

    function [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nz)
    % The function computes the Lagrange-multiplies of the boundary
    % constraints, given an unconstrained field in the frequency domain,
    % evaluated at the omega values which participate in the optimization.
        % Computing the unconstrained field at t=0, by the DCT of the
        % driving field specrum evaluated at t=0:
        fieldt_unc0 = sqrt(2/pi)*sum(fieldw_unc_nz.*integwnz);
        % Computing the unconstrained field at t=T, by the DCT of the
        % driving field specrum evaluated at t=0:
        fieldt_uncT = sqrt(2/pi)*sum(fieldw_unc_nz.*coswT(iEnz).*integwnz);
        lambda0 = (fieldt_unc0*dctfilterE0 - fieldt_uncT*dctfilterET)/deTdctfilterE;
        lambdaT = (fieldt_uncT*dctfilterE0 - fieldt_unc0*dctfilterET)/deTdctfilterE;        
    end

    function fieldw_con = getfieldwcon(fieldw_unc_nz)
    % The function computes a constrained field spectrum, fieldw_con, with 0 boundaries in the
    % time domain, from the unconstrained field spectrum fieldw_unc.
        [lambda0, lambdaT] = fieldw_unc2lambda(fieldw_unc_nz);
        fieldw_con = fieldw_unc_nz - vfilterEnz.*(lambda0 + lambdaT*coswT(iEnz));
    end

    function get_chiihterm(allapsi)
    % The function computes the inhomogeneous term for the chi propagation
    % by the inhomogeneous Schreodinger equation.
        % The filtered evaw, transformed to the time grid and evaluated at
        % the semi-global propagation grid:
        evafil = dctIintgrid(evaw.*vfiltera, T, t_ts(1:(Nt_ts-1)))/dctfactor;
        % The inhomogeneous term at the propagation grid:
        chiihterm = -allapsi(:, allt_lasti:-1:1)*spdiags(evafil(allt_lasti:-1:1), 0, allt_lasti, allt_lasti);     
    end

end