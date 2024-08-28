function [out_sig] = gen_random_pi_sig(Period, NumChannels,...
    NumPeriod, type, band)
% Generate a random signal of NumChannels dimensions and of length Period 
% wich is repeated NumPeriod times. The obtained multidimensional signal is 
% hence periodic and contains Period*NumPeriod samples. 
%
% In: - Period: of the periodic signal to generate (in nb of samples).
%     - NumChannels: dimension of the generated signal. 
%     - NumPeriod: nb of periods of the signal to generate.
%     - type: type of distribution for the random values observed in each
%             period, among the following:
%            'randn': gaussian samples
%            'rand': uniform samples on [-1, 1].
%    - band: (1x2) vector defining the frequency band to consider.
%            Frequencies should be normalized between 0 and 1 so that 1
%            corresponds to the Nyquist frequency (half the sample rate).
%            Default: band =[0 1].
%
% Out: out_sig: (NumChannels x NumPeriod*Period) matrix containing the 
%               generated signals. 
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be


ideal_filter = 1 ;      % otherwise: use a fir filter
plot_filtering = 0 ;    % show the filtered signal (if ideal_filter).

if nargin<5, band = []; end
if nargin<4, type = 'randn' ;end

if isempty(band),band=[0 1];end

Nsamples = NumPeriod*Period ; 

if strcmpi(type,'randn')
    out_sig = randn(NumChannels, Period) ; 
elseif strcmpi(type,'rand')
    % interval considered: [-1, 1]
    out_sig = -1 + 2.*rand(NumChannels, Period) ; 
else
   error(['-- Distribution ', type, ' is unknown']) 
end

out_sig = repmat(out_sig, 1, NumPeriod) ; 
out_sig = out_sig-repmat(mean(out_sig,2), 1, Nsamples) ; 

if ~isequal(band, [0,1])
    if ideal_filter            
        ft_outsig = fft(out_sig,[],2) ; % fft over each row

        % normalized frequencies between 0 and 1, 1 corresponding to
        % fs/2 (nyquist freq)
        freqs=[0:Nsamples-1].*(2/Nsamples) ; 
        freqs(freqs >= 1) = freqs(freqs >= 1) - 2;

        ft_outsig_filt = ft_outsig ; 
        ft_outsig_filt(:,abs(freqs)<band(1) | abs(freqs)>band(2)) = 0 ; 
        outsig_filt = real(ifft(ft_outsig_filt,[],2)) ;

        if plot_filtering
            figure()
            subplot(131)
            plot(1:Nsamples, out_sig(1,:)) ; hold on

            subplot(132)
            plot(freqs, abs(ft_outsig)) ; hold on
            title('Original spectrum')


            subplot(133)
            plot(freqs, abs(ft_outsig_filt)) ; hold on
            title('Filtered spectrum')

            subplot(131)
            plot(1:Nsamples, ifft(ft_outsig,[],2)) ; hold on
            plot(1:Nsamples, outsig_filt) ; hold on
            legend('Original', 'ifft-fft','Filtered')
        end
        out_sig = outsig_filt ; 
    else
        % bandpass-filter the signal
        max_filt_order = 100 ; % cfr PAC_analysis (compare_filtering)
        type_fir1 = 'bandpass' ; 
        if abs(band(1)-0)<10^(-10) 
            type_fir1 = 'low' ;
            band = band(2) ; 
        elseif abs(band(2)-1)<10^(-10)
            type_fir1 = 'high' ; 
            band = band(1) ; 
        end

        filt_order = min(max_filt_order, floor(Nsamples/3)) ;
        b = butter(filt_order,band,type_fir1) ; % or butter (iir)...
        for idx_chan=1:NumChannels
            % Zero-phase digital filtering (a=1)
            out_sig(idx_chan, :) = filtfilt(b, 1, double(out_sig(idx_chan, :))) ;         
        end
    end
end

end

