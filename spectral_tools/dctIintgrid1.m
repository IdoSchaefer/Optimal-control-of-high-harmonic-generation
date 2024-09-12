function dctv = dctIintgrid1(v, interval, int_grid)
% The function computes the dctI of v (see details in the function dctI), 
% in an internal grid within the intervals of the equally spaced 
% Fourier grid.
% The procedure is based on forward time-propagation of the signal in the
% Fourier frequency-domain for all the required time-points in the internal
% grid, and transformation to the time-domain.
% interval:  The length of the whole interval of v in the new space (time space).
% int_grid:  The internal grid inside the interval between equally spaced points.
% It is assumed that v and int_grid are either both column vectors, or
% are both row vectors.
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
    % The normalized internal grid:
    if dim(1) == 1
        normalig = int_grid/dt;
    else
        normalig = int_grid.'/dt;
    end
    % A symmetric extension of v for the representation of v in the Fourier domain:
    v_ext = [v(1:N); v((N - 1):-1:2)];
    M = v_ext*ones(1, Nig);
    % Time-propagation in the Fourier frequency-domain:
    M([1:(N - 1), (N + 1):Nfourier], :) = M([1:(N - 1), (N + 1):Nfourier], :).*exp(-1i*wdt([1:(N - 1), (N + 1):Nfourier])*normalig);
    M(N, :) = M(N, :).*cos(pi*normalig);
    % Transformation to the time-domain:
    fftM = fft(M);
    for gridi = 1:Nig
        dctv(gridi:Nig:(Nfourier*Nig)) = fftM(:, gridi);
    end    
    dctv = dctv(1:((N - 1)*Nig + 1))/sqrt(Nfourier);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end