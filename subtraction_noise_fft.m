function [tf_dataset] = subtraction_noise_fft(tf_dataset, bins_noise)
% Remove the noise from fft data.
%
% Used in TCS2_analyze_SSEP.m
%
% In: - tf_dataset: (n_chan x n_freqs) dataset of fft coefficients
%     - bins_noise: vector indicating the bins to consider to compute the
%                   average noise amplitude in the fft coefficients.
 
[n_chan, n_freqs] = size(tf_dataset) ; 
n_bins = length(bins_noise) ; 

% Average noise amplitude for each entry of tf_dataset
tf_noise = zeros(n_chan, n_freqs) ; 

% Nb of bins considered to compute the average noise for each entry of
% tf_dataset
n_bins_all_tf = n_bins*ones(n_chan, n_freqs) ;

for idx_shift=1:n_bins
    curr_bin = bins_noise(idx_shift) ; 
    
    if curr_bin<0
        tf_noise(:,(1-curr_bin):end) = ...
            tf_noise(:,(1-curr_bin):end) + ...
            tf_dataset(:,1:(end+curr_bin)) ;
        
        n_bins_all_tf(:, 1:(-curr_bin)) = ...
            n_bins_all_tf(:,1:(-curr_bin)) - 1 ; 
        
    else % curr_bin>=0
        tf_noise(:,1:(end-curr_bin)) = ...
            tf_noise(:,1:(end-curr_bin)) + ...
            tf_dataset(:,(1+curr_bin):end) ;
        
        n_bins_all_tf(:,(end-curr_bin+1):end) = ...
            n_bins_all_tf(:,(end-curr_bin+1):end) - 1 ;
    end
    
end
tf_dataset = tf_dataset - tf_noise./n_bins_all_tf ; 

end
