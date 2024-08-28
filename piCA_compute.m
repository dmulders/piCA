function [ W, D] = piCA_compute(dataset, nH, version_piCA, ...
    freqs_vec_to_max, fs)
% Compute un unmixing matrix following one of the criteria of periodic
% component analysis (piCA).

% In: 
%   dataset:   structure with at least the fields:
%               'samples': containing a (C,T) matrix, where C is the
%               number of channels and T the number of time steps of the
%               recordings.  
%               'period': period (in number of samples) of the activity
%               to highlight.
%   nH:     number of harmonics to consider, if version_piCA==3. 2*nH
%           signals will be considered for the reference periodic signal.
%   version_piCA:       in {1,2,3, 4}
%                       1: piCA v1
%                       2: piCA v2
%                       3: CCA approach, cfr Nakanishi 2014.
%                       4: spectral contrast maximization (Sameni). 
%   freqs_vec_to_max:   vector of frequencies (Hz) to enhance, if
%                       version_piCA=4. Should contain positive and
%                       negative frequencies to ensure realness of GEVD.
%   fs:         sampling frequency when version_piCA=4. 
%
% Out:
%   W:      matrix of filters (one filter in each column) found according
%           to the method specified in argument.
%           Columns of W are ranked in decreasing order of periodicity. 
%           ATTENTION: there might be less columns than the input number of
%           channels, since too small eigenvalues are discarded. 
%           Will have to consider the pseudo-inverse to obtain the spatial
%           patterns (pinv(W')).
%   D:      diagonal matrix containing the eigenvalue associated to each
%           column of W.
%
% 
% Dounia Mulders - dounia.mulders@uclouvain.be


FP_weighting = 0 ; % NOT a good idea (at least on real data, for piCA v4)  
% weight the fft coefficient at the harmonics by the corresponding
% frequency to account for the power law decrease in EEG amplitude (piCA v4)

if nargin<4
    if version_piCA==4
        error(['Cannot apply method ',num2str(version_piCA),' without a ',...
            'frequency vector to amplify and sampling rate'])
    end
    freqs_vec_to_max = [] ; 
    fs = 0 ; 
end
if nargin<3
   version_piCA = 1 ;  
end
if nargin<2
    nH = 4 ;
end

samples = dataset.samples ;
[n_chan, nT1, ~] = size(samples) ;

% Fundamental period
TSS = dataset.period ; 

% Standardize
samples = standardize_signals(samples, 0);

ignore_borders = 1 ;    % do not use xcorr by cyclically completing signals

% Periodic Component Analysis (\pi CA)            
switch version_piCA
    case 1        
        % (1) min_w (w' A1_TSS w)/(w' C1 w)
        % Build the matrix A_1(T_{SS})                
        A1_TSS = get_diff_lag_corr_mat(samples, TSS,...
            ignore_borders) ; 

        % Correlation matrix (denominator)
        C1 = get_corr_mat(samples) ;

        [W,D] = eig(A1_TSS, C1) ;
        [W,D] = sort_W_based_on_D(W,D, 'ascend') ;
        
    case 2
        % \pi CA, Sameni 2011
        % (2) max_w (w' C1(T_{SS}) w)/(w' C1 w)                
        C1_TSS = get_xcorr_lag_mat(samples, TSS,...
            ignore_borders) ;  

        % Correlation matrix (denominator)
        C1 = get_corr_mat(samples) ;

        [W,D] = eig(C1_TSS, C1) ;
        [W,D] = sort_W_based_on_D(W,D, 'descend') ;
                
    case 3    
        % (3) CCA
        Y = zeros(2*nH, nT1) ; 
        time_vec = 2*pi/TSS*[0:(nT1-1)] ; 
        idx_row = 1 ; 
        for idxH = 1:nH
            Y(idx_row,:) = sin(time_vec*idxH) ; 
            Y(idx_row+1,:) = cos(time_vec*idxH) ; 
            idx_row = idx_row + 2 ; 
        end

        [W,Wy,D] = canoncorr(samples', Y') ; 
        D = diag(D) ; 
        [W,D] = sort_W_based_on_D(W,D, 'descend') ; 
                
    case 4
        % (4) Spectral contrast maximization                
        C_freq = get_xspectrum_mat(samples, ...
            fs, freqs_vec_to_max, FP_weighting) ;
        C1 = get_corr_mat(samples) ; 

        [W,D] = eig(C_freq, C1) ;
        [W,D] = sort_W_based_on_D(W,D, 'descend') ; % maximize

end
% normalize the spatial filters
W = W./(repmat(sqrt(sum(W.^2, 1)), n_chan,1)) ; 

end


function [sub_tf_samples,freqs_samples] = keep_freq_components(tf_samples, freqs, freqs_vec_to_max)
% From the fft tf_samples (n_chan, n_freqs) whose frequencies are given in
% the vector freqs, return a subsample of this fft with only the
% frequencies from freqs_vec_to_max (negative and positive).

freq_resol = abs(diff(freqs(1:2))) ; 

n_freq_to_max = length(freqs_vec_to_max) ; 
n_chan = size(tf_samples, 1) ; 
sub_tf_samples = zeros(n_chan, n_freq_to_max) ; 
freqs_samples = [] ;  
% negative freq should also be in freqs_vec_to_max!
idx_freq_kept = 0 ;

for idx_freq_max = 1:n_freq_to_max
    curr_freq_max = freqs_vec_to_max(idx_freq_max) ;
    [~, curr_idx] = min(abs( freqs- curr_freq_max)) ;
    freq_in_freqs = freqs(curr_idx) ;
    
    if abs(freq_in_freqs-curr_freq_max)<=freq_resol
        idx_freq_kept = idx_freq_kept + 1 ;
        % keep this frequency, it is in the range covered by freqs
        sub_tf_samples(:, idx_freq_kept) = tf_samples(:, curr_idx) ;
        freqs_samples = [freqs_samples, curr_freq_max] ; 
    end
    
end
sub_tf_samples = sub_tf_samples(:, 1:idx_freq_kept) ; 

end

function xspectrum_matrix = get_xspectrum_mat(data_tensor, ...
    fs, freqs_vec_to_max, FP_weighting)
% Return the cross spectrum matrix computed using the frequencies
% freqs_vec_to_max for the (n_chan, nTimes, n_epochs) dataset data_tensor.
%
% Inputs
%   data_tensor:    (n_chan, nTimes, n_epochs) data set. n_epochs can be 
%                   = 1.
%   fs:     sampling frequency.
%   freqs_vec_to_max:       vector of frequencies at which the cross
%                           spectrum matrix will be computed.
%   FP_weighting:       weight the fft amplitudes by their corresponding
%                       frequencies to accoutn for the power law decrease
%                       of EEG sepctra. 
%

if nargin<4
   FP_weighting = 0 ;  
end
[n_chan,~,n_epochs] = size(data_tensor) ; 
n_freq_max = 0 ; % total, for all the epochs
xspectrum_matrix = zeros(n_chan, n_chan) ; 

% normalization factor used when FP_weighting
norm_all_freq = max(1,min(abs(freqs_vec_to_max))) ;

for idx_epoch=1:n_epochs
    samples_tmp = squeeze(data_tensor(:,:,idx_epoch)) ;
    [~, freqs, tf_samples] = get_spectrum(samples_tmp', fs, 0, ...
        1) ;
    tf_samples = tf_samples' ; 
    [tf_samples, freqs_samples] = keep_freq_components(tf_samples, freqs, ...
        freqs_vec_to_max) ; 
    curr_n_freq_max = size(tf_samples, 2) ; 
    n_freq_max = n_freq_max + curr_n_freq_max ;
    
    if FP_weighting
       freqs_samples = abs(freqs_samples./norm_all_freq) ;  
       tf_samples = tf_samples.*repmat(reshape(freqs_samples, 1, ...
           curr_n_freq_max), n_chan,1) ; 
    end    
    
    % average the cov matrices of all the trials
        
    % real: because otherwise it might remain small complex
    % parts (10^(-18) imag parts)
    xspectrum_matrix = xspectrum_matrix + real(((tf_samples*tf_samples')...
        ./curr_n_freq_max)./n_epochs) ;
end

end

function corr_matrix = get_corr_mat(data_tensor)
% Return the multidimensionnal correlation matrix for the signals
% data_tensor.
%
% Inputs
%   data_tensor:    (n_chan, nTimes, n_epochs) data set. n_epochs can be 
%                   = 1. Each signal is supposed to be centered.
%

[n_chan, nT, n_epochs] = size(data_tensor) ; 

% average the corr matrix of all trials
corr_matrix = zeros(n_chan, n_chan) ;
for idx_epoch=1:n_epochs
    samples_tmp = squeeze(data_tensor(:,:,idx_epoch)) ;
    % average the cov matrices of all the trials
    corr_matrix = corr_matrix + ((samples_tmp*samples_tmp')./nT)./n_epochs ;                            
end    

end

function diff_n_lag_corr = get_diff_lag_corr_mat(data_tensor, n_lag,...
    ignore_borders)
% Return the correlation matrix of the difference of the input signal with
% itself delayed from n_lag samples.
%
% Inputs
%   data_tensor:    (n_chan, n_times, n_epochs) data set. 
%   n_lag:          number of samples used for the lag.
%

% Build the matrix A(n_lag)

if nargin<3
   ignore_borders = 1 ;  
end

[n_chan, nT, n_epochs] = size(data_tensor) ; 

   
if ignore_borders                        
    diff_n_lag = data_tensor(:,1:(end-n_lag), :) - ...
        data_tensor(:, (1+n_lag):end,:) ;
    % center the signals
    diff_n_lag = standardize_signals(diff_n_lag, 1) ;

    nT_used = nT - n_lag ;
else
    diff_n_lag = data_tensor - ...
        data_tensor(:, [(1+n_lag):end, 1:n_lag],:) ;
    nT_used = nT ; 
end

diff_n_lag_corr = zeros(n_chan, n_chan) ; 
for idx_epoch=1:n_epochs
    diff_tmp = squeeze(diff_n_lag(:,:,idx_epoch)) ; 
    diff_n_lag_corr = diff_n_lag_corr + ((diff_tmp*diff_tmp')./nT_used)./n_epochs ; 
end


end

function corr_n_lag = get_xcorr_lag_mat(data_tensor, n_lag,...
    ignore_borders)
% Compute the cross correlation matrix between the input signal and its
% shifted version by n_lag.

if nargin<3
   ignore_borders = 1 ;  
end

[n_chan, nT, n_epochs] = size(data_tensor) ; 

corr_n_lag = zeros(n_chan, n_chan) ; 
for idx_epoch=1:n_epochs
    samples_tmp = squeeze(data_tensor(:,:,idx_epoch)) ; 
    if ignore_borders
        % center all the used signals
        samples_tmp_first = standardize_signals(...
            samples_tmp(:,1:(end-n_lag)), 1) ;
        samples_tmp_second = standardize_signals(...
            samples_tmp(:, (1+n_lag):end), 1) ;

        corr_n_lag = corr_n_lag + (samples_tmp_first*...
            samples_tmp_second')./((nT-n_lag)*n_epochs) ;

    else
        corr_n_lag = corr_n_lag + (samples_tmp*...
            (samples_tmp(:, [(1+n_lag):end, 1:n_lag]))')./(nT*n_epochs) ; 
    end
end
    

corr_n_lag = (corr_n_lag + corr_n_lag')./2 ; % ensure symmetry

end
