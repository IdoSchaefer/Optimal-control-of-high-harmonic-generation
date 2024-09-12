load coulomb_optV240
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt15, fieldw15, psi15, evat15, evaw15, evmut15, evmiuw15, relE15, conv15, niter15, mallniterc15, J115, maxgrad15, alpha15, invHess15] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
% Reset of the Hessian approximation after the optimization gets stuck:
[fieldt15a, fieldw15a, psi15a, evat15a, evaw15a, evmut15a, evmuw15a, relE15a, conv15a, niter15a, mallniterc15a, J115a, maxgrad15a, alpha15a, invHess15a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw15, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.9).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
