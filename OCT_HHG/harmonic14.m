load coulomb_optV240
options = optionsOCqn(1e-4, 1e3);
options.f_max_alpha = get_f_max_alphaOCf(0.15, 0.2, 1e3, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)));
[fieldt14, fieldw14, psi14, evat14, evaw14, evmut14, evmuw14, relE14, conv14, niter14, mallniterc14, J114, maxgrad14, alpha14, invHess14] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, @(w) 5*exp(-(w-0.06).^2/(2*0.01^2)).*sin((w-0.06)*pi/0.015), @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.84).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
% Reset of the Hessian approximation after the optimization gets stuck:
[fieldt14a, fieldw14a, psi14a, evat14a, evaw14a, evmut14a, evmuw14a, relE14a, conv14a, niter14a, mallniterc14a, J114a, maxgrad14a, alpha14a, invHess14a] = OCfx_qn(fi0240, Vabs240, 1, [-240 240], xabs240, -x240./(1 + x240.^2).^(3/2), @(x) 0.01*0.5*(tanh(50*(x-0.9)) - tanh(5)), @(x) 0.01*0.5*50*sech(50*(x-0.9))^2, fieldw14, @(w) 5e5*exp(-(w-0.06).^2/(2*0.01^2)), @(w) exp(-(w-0.84).^2/(2*0.01^2)), options, 1e3, 0.2, 7, 7, 1e-4, 1e3);
