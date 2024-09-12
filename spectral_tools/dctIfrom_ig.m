function dctv = dctIfrom_ig(v, interval, int_grid, weights)
% The inverse of the function dctIintgrid;
% The function computes the DCT-I of a time signal v, which is
% sampled in a periodic grid structure with an internal grid in each
% time-step. The time-step length is uniform.
% interval:  The whole interval of v.
% int_grid:  The internal grid inside the interval between equally spaced
% points, in the time space.
% weights: The integration weights of the points in the internal grid.
%%% Note: This program can be improved---there are operations which can be
%%% performed without a loop, and there are also unnecessary operations (see
%%% comment below). In addition, the procedure can be simplified and made more efficient if the
%%% internal grid has a symmetric structure around the central point, which
%%% is the common case (e.g. Chebyshev sampling). There is an improved
%%% version for a symmetric internal structure in the procedure
%%% dctIfrom_ig_sym1.m, with more detailed comments.
    Nig = length(int_grid);
    dim = size(v);
    if dim(1) == 1
        Np = dim(2);
        v = v.';
    else
        Np = dim(1);
    end
    N = (Np - 1)/Nig + 1;
    Nfourier = 2*(N - 1);
    dctv = weights(1)*ifft([v(1:Nig:Np); v((Np - Nig):-Nig:(Nig + 1))]);
    wdt = (0:2*pi/Nfourier:(2*pi*(1 - 1/Nfourier))).';
    % Note that the points of wdt after N are irrelevant, see the computation of dctv.
    % I mistakenly didn't substract 2*pi from the last Nfourier/2 terms,
    % but this doesn't matter for the result.
    dt = interval/(N - 1);
    normalig = int_grid/dt;
    vgridi = zeros(Nfourier, 1);
    for gridi = 2:Nig
        vgridi(1:(N - 1)) = v(gridi:Nig:Np);
        dctv = dctv + weights(gridi)*(exp(1i*wdt*normalig(gridi)).*ifft(vgridi) + ...
            exp(-1i*wdt*normalig(gridi)).*ifft([vgridi(1); vgridi(Nfourier:-1:2)]));
    end
    dctv = dctv(1:N)*sqrt(Nfourier);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end