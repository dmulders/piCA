function [xtf, freqs, xtf_complex] = get_spectrum(x, f_s, only_pos, divide_by_N,...
    padd_fft,N_padd)
% Returns the spectrum of the input signal, together with the associated
% frequencies.
% if x is a matrix, the fft is computed along each column
%

% used in test_EEG_model

if nargin<5
    padd_fft = 0 ; 
    N_padd = NaN ; 
end
if nargin<4
    % normalization
   divide_by_N = 0 ;  
end
if nargin<3
    only_pos = 0 ; % only consider positive frequencies (eg if real signal)
end
 
if isrow(x)
    % ensure that x is a column vector
    x = x' ; %reshape(x,length(x),1) ;    
end
Sx = size(x) ;
Nfft = Sx(1) ; 
if padd_fft
    xtf = fft(x,N_padd) ;
else
    xtf = fft(x) ;
end
N_freqs = size(xtf,1) ; 
% fft of each column; corresponding freqs: [0,fs]
% after fftshift: freqs in [-fs/2,fs/2]
if nargout>2
   xtf_complex = xtf ; 
end

for idx_col=1:Sx(2)
    xtf(:,idx_col) = abs(fftshift(xtf(:,idx_col))) ; 
    if nargout>2
       xtf_complex(:,idx_col) = fftshift(xtf_complex(:,idx_col)) ; 
    end
end
if divide_by_N
    xtf = xtf./Nfft ; 
    if nargout>2
       xtf_complex = xtf_complex./Nfft ; 
    end
end

freqs=[0:N_freqs-1].*(f_s/N_freqs);
freqs(freqs >= f_s/2) = freqs(freqs >= f_s/2) - f_s;
freqs=fftshift(freqs);
% freqs = freqs(1:fix(Nfft/2)) ; 

if only_pos    
    idx_neg_freq = freqs < 0 ; 
    freqs(idx_neg_freq) = [] ;
    xtf(idx_neg_freq,:) = [] ;
    if nargout>2
       xtf_complex(idx_neg_freq,:) = [] ; 
    end
    idx_double_freq = freqs~=0 & freqs~=f_s/2 ; 
    xtf(idx_double_freq',:) = 2*xtf(idx_double_freq',:) ; 
    % double the power of all frequencies except the Nyquist one    
end

if isrow(x)
    % ensure that x is a column vector
    xtf = xtf' ; %reshape(xtf,Sx_orig) ;    
    if nargout>2
       xtf_complex = xtf_complex' ; 
    end
    
end


end