function [simu_data] = generate_mixed_sig(NF, stds_or_N, SNR, nChan, ...
    fs, fZ, nPeriod, nHarmonics, inv_freq_power, opt_std, rng_opt)
% Stimulate data according to a linear mixing model of periodic sources. 
%
% Model:
% .-------------------.
% | X_k = A*Z_k + n_k |
% .-------------------.
%   k: indicates the time sample
%   n_k: noise with inverse frequency power inv_freq_power and with std 
%          defined to reach the input SNR
%
%
% All the "sources" (or factors) in the model are periodic of period 1/fZ.
% The multimodal sources have the same time courses, but they are delayed 
% of a time lag tau for the 1st class. This lag can account for the 
% differences in conduction velocities of the afferent fibers.
%
% In: - NF: nb of hidden factors
%     - stds_or_N: std of signals A_s*Z^{s(i)} or n^i, depending on
%                  opt_std. 
% 
%     - SNR: targeted signal-to-noise ratio, in dB, from 
%             which the variance of the noise (resp the signal) is fixed if
%             opt_std = 1 (resp. = 2).
%     - nChan: number of recording channels
%     - fs: sampling frequency
%     - fZ: frequency of the factors Z (ie SS frequency if we work with
%            steady-state signals).
%     - nPeriod: number of periods (of length 1/f_Z seconds) to generate.
%     - nHarmonics: number of harmonics to keep in the factors Z
%     - inv_freq_power: inverse frequency power of the noise 
%                       (1/f^inv_freq_power spectral density)
%                       Default: 1, pink noise.
%     - opt_std: in {1,2}.
%                 1: stds_or_N is the standard deviation of A_s*Z^{s(i)}. 
%                 2: stds_or_N is the standard deviation of the noise n^i.
%     - rng_opt: seed for the random number generator. Can be 'default',
%                'shuffle' or any integer >= 0. Is not set by default.
%
% Out: - simu_data: structure with fields samples and mixingMat.
%           samples: (channels x timesteps) matrix containing the EEG 
%                     signals.
%           mixingMat: mixing matrix used
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be
%
%%%%%%%%%%%%%%%%%%%%%%%%%
lowpass_noise = 0 ;     % low-pass filter the pink noise before fixing its 
                        % std to reach the targeted SNR.
%%%%%%%%%%%%%%%%%%%%%%%%%
if NF<1
   error('There should be at least one hidden factor.') 
end
if nargin>10
    % if 'Default': reset the random number generator used by RAND, RANDI, 
    % and randn to its default startup settings, so that randn produces 
    % the same random numbers as if you restarted MATLAB.
    rng(rng_opt) ; 
end
if nargin<10
    opt_std = 1 ; 
end
if nargin<9
   inv_freq_power = 1 ; % consider a pink noise (1/f noise) 
end

% 10*log10(var(A_s*Z^{s(i)})/var(n^i)) = SNR 
% --> var(n^i) = 10^(-SNR/10)*var(A_s*Z^{s(i)})
% --> var([A_s*Z^{s(i)}](1,:)) = 10^(SNR/10)*var(n^i)
if opt_std==1  
    stds = stds_or_N ;
    stdN = stds*10^(-SNR/20) ;    
elseif opt_std==2
    stdN = stds_or_N ;
    stds = stdN*10^(SNR/20) ;           
else
    error('--- opt_std can only be 1 or 2')  
end

% *** Generate mixing matrix: uniform entries in [-1, 1] 
% (in  Wu 2009: in [0,1]).
mixingMat = rand(nChan, NF) ;

% normalize so that norm of each column = 1
mixingMat = mixingMat./repmat(sqrt(sum(mixingMat.^2,1)), nChan, 1) ; 


fLastHarmonic = fZ*nHarmonics ;             % Keep nHarmonics harmonics.
band_hz = [0, min(fs/2, fLastHarmonic)] ; 
% normalized: 1 corresponds to the Nyquist frequency = fs/2
band_norm = band_hz./(fs/2) ; 
Period = round(fs/fZ) ; % in nb of samples
nSamples = Period*nPeriod ; 

% *** Generate the factors

% std for each factor
if NF==1
    stds_vec_Zs1 = 1 ; 
else
    stds_vec_Zs1 = [(NF-1):-1:0] ; 
end

Zs1 = gen_random_pi_sig(Period, NF,...
    nPeriod, 'randn', band_norm) ; % one channel per row
Zs1 = Zs1.*repmat(stds_vec_Zs1', 1, nSamples)./repmat(std(Zs1,[],2), 1,nSamples) ; 

% Impose the standard deviation of the "s" factors
tmp_mixed_Zs1 = mixingMat*Zs1 ;
std_Zs1 = std(tmp_mixed_Zs1(:)) ; 


Zs1 = (Zs1.*stds)./std_Zs1 ; 

% *** Factor model for the mixtures
% --> specific factors (always at least 1)
mixings = mixingMat*Zs1 ;% mixingMat*[Zu1; Zs1] ; 

% *** Add the pink noise
cn = dsp.ColoredNoise('InverseFrequencyPower',inv_freq_power,...
    'NumChannels',nChan, 'SamplesPerFrame',nSamples);
noise = (step(cn))' ; % nSamples x nChan

maxFreqNoise = fLastHarmonic*2 ;
% lowpass filter the pink noise below this frequency to concentrate its
% power on lower frequencies
if lowpass_noise && maxFreqNoise<fs/2
    % Lowpass filtering 
    % Attention: EEGLAB has its own "firls" function, which creates
    % conflicts with the original matlab one...
    max_filt_order = 100 ;
    filt_order = min(max_filt_order, floor(nSamples/3)) ;
    b = fir1(filt_order,maxFreqNoise/(fs/2),'low') ;
    for idx_chan=1:nChan
        % Zero-phase digital filtering (a=1)
        noise(idx_chan, :) = filtfilt(b, 1, double(noise(idx_chan,:)));         
    end
end

% impose the standard deviation
noise = noise.*stdN./std(noise(:)) ;
samples = mixings + noise ; 

% Pink noise + mixture
samples = standardize_signals(samples,0) ; % standardize_signals(samples_cls2, 0);

simu_data = struct('samples', samples, 'mixingMat', mixingMat, ...
    'first_source', Zs1(1,:)) ; 


end

