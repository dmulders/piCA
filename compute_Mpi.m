function [Mpi_signal] = compute_Mpi(in_signal, fs, f_pi, n_harm_save, bins_noise)
% Compute the periodicity measure M_{\pi} for the input signal in_signal.
% Reflect global and local properties of spectrum:
% - globally, peaks (freq content) at the harmonics should be present
% - locally, these peaks should be larger than the neighboring amplitudes
%            at the other frequencies.
%
% Inputs:
%   - in_signal: signal for which the periodicity measure has to be
%               computed. Has to be (n_chan, n_times).
%   - fs: sampling frequency (Hz).
%   - f_pi: fundamental frequency of interest (Hz).
%   - n_harm_save: number of harmonics, including the fundamental frequency
%                   that has to be accounted for in the computation of Mpi.
%
% Outputs:
%   - Mpi_signal: periodicity measure for each channel.
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
divide_tf_by_N = 1 ;            %
if nargin<5
    bins_noise = [-5:-2,2:5] ;  % bins (wrto concerned value) of the 
                                % frequency spectrum to consider, if 
                                % subtract_noise_fft    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscolumn(in_signal)
    in_signal = in_signal' ; 
end

n_times = size(in_signal, 2) ; 

% Compute the fourier transform:
[tf_dataset, freqs] = get_spectrum(in_signal', fs, 1, ...
    divide_tf_by_N) ;
tf_dataset = tf_dataset' ; % (n_chan, n_freqs)

norm_factor = sum(tf_dataset, 2) ; 
tf_dataset = 100.*tf_dataset./(repmat(norm_factor, ...
    1, size(tf_dataset, 2))) ; 

if length(bins_noise)>0
    tf_dataset_NS = subtraction_noise_fft(tf_dataset, bins_noise) ;
else
    tf_dataset_NS = tf_dataset ; % without any noise-subtraction
end
% Denominator of Mpi:
ampl_signal = sum(tf_dataset, 2) ; % all equal 100 if tf already normalized above...
% Numerator of Mpi:
all_nH_values = zeros(size(ampl_signal)) ; % (nChan, 1)

% Sum of all harmonics up to max_freq_harm
n_poss_harm = floor(freqs(end)/f_pi) ; 
n_harm_used = min(n_harm_save, n_poss_harm) ; 

for idx_harm=1:n_harm_used                        
    curr_freq = idx_harm*f_pi ;                         
    [~,idx_freq] = min(abs(freqs-curr_freq));

    all_nH_values = all_nH_values + tf_dataset_NS(:,idx_freq) ; 
end 


% One value per signal component (if multidim)
Mpi_signal = 100*all_nH_values./ampl_signal ; 

end

