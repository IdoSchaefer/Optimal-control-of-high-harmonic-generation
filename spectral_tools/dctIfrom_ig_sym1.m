function dctv = dctIfrom_ig_sym1(v, interval, int_grid, weights)
% The inverse of the function dctIintgrid; applies for a symmetric internal
% grid structure.
% The function computes the DCT-I of a time signal v, which is
% sampled in a periodic grid structure with a symmetric internal grid in each
% time-step. The time-step length is uniform.
% The integration is based on a weighted sum of the points in the internal
% grid. The FFT procedure is employed for efficiency, by the decomposition of the time-signal
% into several signals with equally-spaced sampling, for each of the internal points. Each
% signal is transformed to the Fourier frequency-domain separately by an iFFT, and
% backward propagated in time to the beginning of the time-interval. The
% different resulting frequency signals are used to obtain the result by a
% weighted sum according to the integration weights for each frequency
% component.
% This practice can be shown to hold by a detailed mathematical analysis (still in my notebook, 28.1.20).
% interval:  The whole interval of v.
% int_grid:  The internal grid inside the interval between equally spaced
% points, in the time space.
% weights: The integration weights of the points in the internal grid.
% It is assumed that v, int_grid and weights are either all column vectors, or
% are all row vectors.
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
    % The symmetric extension of v, for representation in the Fourier
    % domain:
    v_ext = [v(1:Np); v((Np - 1):-1:2)];
    Next = 2*Np - 2;
    Mv_ext = zeros(Nfourier, Nig);
    % Decomposition of the signal into several signal with equally spaced
    % sampling:
    for gridi = 1:Nig
        Mv_ext(:, gridi) = v_ext(gridi:Nig:Next);
    end
    % Transformation of the signals to the Fourier frequency-domain:
    ifft_v_ext = ifft(Mv_ext);
    wdt = (0:2*pi/Nfourier:pi).';
    dt = interval/(N - 1);
    normalig = int_grid/dt;
    % Backward time-propagation in the frequency domain, and a weighted sum
    % for each frequency component:
    if dim(1) == 1
        dctv = sum((ifft_v_ext(1:N, :).*exp(1i*wdt*normalig))*weights.', 2);
    else
        dctv = sum((ifft_v_ext(1:N, :).*exp(1i*wdt*normalig.'))*weights, 2);
    end
    dctv = dctv*sqrt(Nfourier);
    if dim(1) == 1
        dctv = dctv.';
    end
    if isreal(v)
        dctv = real(dctv);
    elseif isreal(1i*v)
        dctv = 1i*imag(dctv);
    end
end