function max_alpha = alpha_maxOCf(fieldwnz0, direction, max_fieldt, Nt, iEnz, T)
% Intended for problems with restriction of the frequency domain of the control field 
% (the OCf_qn programs). 
% The function returns the maximal allowed alpha value for the line-search
% procedure line_search1, for a given maximal allowed value of the 
% magnitude of the control field in the time domain. It is possible to set a uniform
% value for all the transformed x (the infinity norm of the dctI of x), or to specify the
% value for each of the components of the transformed x.
% fieldwnz0: The initial point of the search, i.e. the solution from the previous
% iteration
% direction: The direction of search
% max_fieldt: The value of the maximal allowed magnitude of the field in 
% the time domain; if all the time-points are restricted to a uniform maximal
% magnitude, max_fieldt may be a positive scalar. If there are different maximal
% magnitudes for the different time-points, max_dctx is a vector of
% dimension (Nt + 1).
% Nt: The number of time-points
% iEnz: The indices of non-zero field in the frequency domain (see  the
% OCf_qn programs).
% T: The final time
    direction_with0 = zeros(Nt + 1, 1);
    direction_with0(iEnz) = direction;
    dct_direction = dctI(direction_with0)*sqrt(Nt*pi)/T;
    fieldw0 = zeros(Nt + 1, 1);
    fieldw0(iEnz) = fieldwnz0;
    dct_x0 = dctI(fieldw0)*sqrt(Nt*pi)/T;
    max_alpha = min((max_fieldt - sign(dct_direction).*dct_x0)./abs(dct_direction));
end