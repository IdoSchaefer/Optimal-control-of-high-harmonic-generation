function dctv = dctIintgrid(v, interval, int_grid)
% The function computes the dctI of v (see details in the function dctI), 
% in an internal grid within the intervals of the equally spaced 
% Fourier grid.
% interval:  The length of the whole interval of v in the new space (time space).
% int_grid:  The internal grid inside the interval between equally spaced points.
%%% An improved version with more detailed comments can be found in the
%%% procedure dctIintgrid.m.
    Nig = length(int_grid);
    dim = size(v);
    if dim(1) == 1
        N = dim(2);
        v = v.';
    else
        N = dim(1);
    end
    dt = interval/(N - 1);
    Nfourier = 2*(N - 1);
    wdt = (0:2*pi/Nfourier:(2*pi*(1 - 1/Nfourier))).';
    wdt((N + 1):Nfourier) = wdt((N + 1):Nfourier) - 2*pi;
    dctv = zeros(Nfourier*Nig, 1);
    normalig = int_grid/dt;
    for gridi = 1:Nig
        dctv(gridi:Nig:(Nfourier*Nig)) = 1/sqrt(Nfourier)*fft([exp(-1i*wdt(1:(N - 1))*normalig(gridi)).*v(1:(N - 1)); ...
            cos(pi*normalig(gridi))*v(N); exp(-1i*wdt((N + 1):Nfourier)*normalig(gridi)).*v((N - 1):-1:2)]);
    end    
    dctv = dctv(1:((N - 1)*Nig + 1));
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end