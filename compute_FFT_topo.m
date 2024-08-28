function topo_fft = compute_FFT_topo(signals_in, fs, f_in)
% Gives the FFT values at f_in for the multidimensionnal input signal
% signals_in. 
%
% Inputs
%   - signals_in:   (nTimes, nChan) signals.
%   - fs:       sampling frequency (Hz).
%   - f_in:     frequency of interest for which topography has to be
%               computed.
%
% Output
%   - topo_fft:     FFT at f_in on all the channels (column vector). 
%

[abs_tf_samples, freqs] = get_spectrum(signals_in, fs, 1, ...
    1) ;
% get_spectrum(x, f_s, only_pos, divide_by_N)
abs_tf_samples = abs_tf_samples' ;
[~,idx_f_in] = min(abs(freqs-f_in)) ; 

topo_fft = abs_tf_samples(:,idx_f_in) ;

end

