function [xalpha, falpha, fgrad_xalpha, sol_achieved, nfevals, alpha] = line_search1(fDf, x0, f0, fgrad0, direction,...
    alpha1, ro, sigma, tau1, tau2, tau3, minf, max_alpha, alpha_factor, Lsection_factor)
% The function implements the line-search procedure described in "Practical
% Methods of Optimization" by Fletcher, Sec. 2.6, for available gradient
% data.
% The difference from line_search is in the definition of the input
% function handle.
%%% Input: %%%
% fDf: a function handle which returns the minimized function and its
% gradient as two different output arguments in the following way:
% [function, gradient] = fDf(x)
% x0: The starting point of the line search
% f0: The f value in the starting point of the line search
% fgrad0: The gradient value in the starting point of the line search
% direction: a vector which represents the direction of search
% alpha1, ro, sigma, tau1, tau2, tau3: Parameters which are defined in the
% reference above
% minf: Defined as \bar{f} in the reference above; if irrelevant, substitute -inf.
% max_alpha: The maximal allowed alpha value (alpha is defined in the
% reference above.)
% alpha_factor: The maximal allowed ratio of the final alpha resulting in
% the bracketing phase, and the initial alpha, represented by alpha1.
% Lsection_factor: The minimal allowed ratio of the final length of the
% sectioning interval and its initial length.
%%% Output: %%%
% xalpha: The x value which is the result of the line search
% falpha: The function value at xalpha
% fgrad_xalpha: The gradient value at xalpha
% sol_achieved: a logical variable which indicates if an acceptable
% solution has been achieved; relevant only if minf is supplied (not
% substituted with -inf).
% nfevals: The number of function and gradient evaluations
% alpha: The resulting alpha in the process
    sol_achieved = false;
    accepted_point = false;
    accepted_interval = false;
    %%% Bracketing phase: %%%
    alpha = alpha1;
    df0_dalpha = fgrad0.'*direction;
    mu = min([abs((minf - f0)/(ro*df0_dalpha)), max_alpha]);
    % The absolute value is taken to account for the case of positive
    % df0_dalpha, which may result from errors in the gradient computation
    % or roundoff errors.
    alpha_last = 0;
    falpha_last = f0;
    Dfalpha_last = df0_dalpha;
    niter_bracket = 0;
    while ~accepted_interval && alpha/alpha1<=alpha_factor && alpha_last<mu && ~isnan(falpha_last)
        xalpha = x0 + alpha*direction;
        [falpha, fgrad_xalpha] = fDf(xalpha);
        Dfalpha = fgrad_xalpha.'*direction;
        if falpha<=minf
            sol_achieved = true;
            accepted_point = true;
            accepted_interval = true;
        elseif falpha>f0 + ro*alpha*df0_dalpha || falpha>=falpha_last
            a = alpha_last;
            b = alpha;
            fa = falpha_last;
            Dfa = Dfalpha_last;
            fb = falpha;
            Dfb = Dfalpha;
            accepted_interval = true;
        elseif abs(Dfalpha)<=-sigma*df0_dalpha
            accepted_point = true;
            accepted_interval = true;
        elseif Dfalpha>=0
            a = alpha;
            b = alpha_last;
            fa = falpha;
            Dfa = Dfalpha;
            fb = falpha_last;
            Dfb = Dfalpha_last;
            accepted_interval = true;
        else
            if 2*alpha - alpha_last<mu
                min_interval = [2*alpha - alpha_last, min(mu, alpha + tau1*(alpha - alpha_last))];
                alpha_new = min_cube(alpha_last, alpha, falpha_last, Dfalpha_last, falpha, Dfalpha, min_interval);
                alpha_last = alpha;
                alpha = alpha_new;
            else
                alpha_last = alpha;
                alpha = mu;
            end
            falpha_last = falpha;
            Dfalpha_last = Dfalpha;
        end
        niter_bracket = niter_bracket + 1;
    end
    if ~accepted_interval
        fprintf('\nWarning: Bracketing phase failed.\n')
        if alpha_last==mu
            fprintf('Search interval was limited by the user preferences.\n')
        end
        if isnan(falpha_last)
            fprintf('The function value is NaN or Inf.\n')
        end
        nfevals = niter_bracket;
        return
    end
    %%% Sectioning phase: %%%
    niter_section = 0;
    if ~accepted_point        
        Lsection = abs(b - a);
        Lsection1 = Lsection;
    end
    while ~accepted_point && 2*(b - a)*Dfa/(abs(fa) + abs(fb))<-10*eps && Lsection/Lsection1>=Lsection_factor
    % The condition 2*(b - a)*Dfa/(abs(fa) + abs(fb))<-10*eps indicates
    % that the estimated relative improvement in the solution in the search is close to
    % the machine precision.
    % (b - a)*Dfa should be <0; its magnitude is an estimation of the upper limit of
    % the magnitude of the absolute improvement in the search.
        min_interval = [a + tau2*(b - a), b - tau3*(b - a)];
        alpha = min_cube(a, b, fa, Dfa, fb, Dfb, min_interval);
        xalpha = x0 + alpha*direction;
        [falpha, fgrad_xalpha] = fDf(xalpha);
        Dfalpha = fgrad_xalpha.'*direction;
        if falpha>f0 + ro*alpha*df0_dalpha || falpha>=fa
            b = alpha;
            fb = falpha;
            Dfb = Dfalpha;
        elseif abs(Dfalpha)<=-sigma*df0_dalpha
            accepted_point = true;
        else
            if (b - a)*Dfalpha>=0
                b = a;
                fb = fa;
                Dfb = Dfa;
            end
            a = alpha;
            fa = falpha;
            Dfa = Dfalpha;
        end
        niter_section = niter_section + 1;
        Lsection = abs(b - a);
    end
    if ~accepted_point
        if 2*(b - a)*Dfa/(abs(fa) + abs(fb))>=-10*eps
            fprintf('\nWarning: Line search failed. The solution is likely to be close to the minimum. An error in the gradient function is also possible.\n')
        else
            fprintf('\nWarning: Sectioning phase failed.\n')
        end
    end
    nfevals = niter_bracket + niter_section;
    % Activate this line in order to check if the search behaves normally:
    %[niter_bracket, niter_section, nfevals]
end

%%%%%% Sub-functions %%%%%%%%%

function minimizer = min_cube(a, b, fa, Dfa, fb, Dfb, interval)
%   The function constructs a cubic interpolation polynomial based on the
%   function data and derivative data in a and b, and returns the minimal
%   value of the interpolation polynomial in the interval indicated by the
%   variable interval.
%   For details see Fletcher (page 37 in the second edition of the reference above).
    L = b - a;
    Dz_fa = Dfa*L;
    Dz_fb = Dfb*L;
    Delta = fb - fa;
    c0 = fa;
    c1 = Dz_fa;
    c2 = 3*Delta - 2*Dz_fa - Dz_fb;
    c3 = Dz_fa + Dz_fb - 2*Delta;
    % The discriminant expression, which appears in the solution to
    % P'(z)=0 (where P(z) is the interpolation polynomial):
    sqrt_expression = c2^2 - 3*c1*c3;
    interval_z = (interval - a)/L;
    % Analysis of the sign of P''(z) shows that the positive root yields a
    % minimum, and the negative yields a maximum. Thus, only the positive root
    % is a candidate for the minimal point.
    if sqrt_expression>0
        zmin = (-c2 + sqrt(sqrt_expression))/(3*c3);
        % According to this specific line-search algorithm, the sign of L
        % is identical to the sign of interval(2)-interval(1), which
        % ensures that interval_z is ordered in increasing values.
        if zmin>interval_z(1) && zmin<interval_z(2)
            % if alpha_min is included in interval:            
            candidates = [interval, a + zmin*L];
            % since alpha_min = a + zmin*L;
            candidates_z = [interval_z, zmin];
        else
            candidates = interval;
            candidates_z = interval_z;
        end
    else
        candidates = interval;
        candidates_z = interval_z;        
    end
    pol_candidates = cube_pol(candidates_z, c0, c1, c2, c3);
    % Actually, this can be made more efficient, since it isn't always necessary
    % to compute the polynomial at the two boundaries. However, it is
    % assumed that the computational time is negligible, and this is
    % much easier to write than a complex system of conditions.
    [~, imin] = min(pol_candidates);
    minimizer = candidates(imin);
end

function result = cube_pol(x, c0, c1, c2, c3)
    result = c0 + c1*x + c2*x.^2 + c3*x.^3;
end
