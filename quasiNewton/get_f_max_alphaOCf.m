function f_max_alphaOCf = get_f_max_alphaOCf(max_fieldt, dt, T, filterE)
% The function returns a function handle f_max_alphaOCf which fits the
% field f_max_alpha in the options structure of the program quasiNewton.m,
% for restriction of the maximal allowed magnitude of the field in the
% time-domain in the OCf_qn programs.
% Input:
% max_fieldt: The value of the maximal allowed magnitude of the field in 
% the time domain; if all the time-points are restricted to a uniform maximal
% magnitude, max_fieldt may be a positive scalar. If there are different maximal
% magnitudes for the different time-points, max_dctx is a vector of
% dimension (Nt + 1).
% dt: The time-step length
% T: The final time
% filterE: As in the OCf_qn programs
    Nt = T/dt;
    w = (0:pi/T:pi/dt).';
    vfilterE = filterE(w);
    iEnz = find(vfilterE>=eps*max(vfilterE));
    f_max_alphaOCf = @(fieldwnz0, direction) alpha_maxOCf(fieldwnz0, direction, max_fieldt, Nt, iEnz, T);
end