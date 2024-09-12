load coulomb_optV240
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt13, fieldw13, psi13, evat13, evaw13, evmut13, evmuw13, relE13, conv13, niter13, mallniterc13, J113, maxgrad13, alpha13, invHess13] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
% Reset of the Hessian approximation after the optimization gets stuck:
[fieldt13a, fieldw13a, psi13a, evat13a, evaw13a, evmut13a, evmuw13a, relE13a, conv13a, niter13a, mallniterc13a, J113a, maxgrad13a, alpha13a, invHess13a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw13, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.78).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
