function [ ] = test_EEG_model( )
% Function to test some models generating synthetic EEG data.
%
% For instance, the following components/features are considered:
% --> pink noise (PSD ~1/f)
% --> periodic random signals (can be used as factors in a factor analysis
%                              model to simulate EEG).
% --> mixture of constrained factor analysis models, one model for each
%       class (eg class "nociceptive" vs "non-nociceptive") 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all ;                %
%rng('default');            %
do_plots = 0 ;              % plots of the pink noise and random period sig
f_s = 500 ;                 %
f_SS = 2 ;                  %
std_pi_sig = 5 ;            % standard deviation for the periodic signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================================= %
% ----------> Generating a pink noise 
% ======================================================================= %
% to model the additive noise
NumChannels = 2;
NumSamples = 1000;
std_noise = 6 ; 
inv_freq_power = 0 ; 
% 0: white noise; 1:pink noise; 2: brownian noise
cn = dsp.ColoredNoise('InverseFrequencyPower',inv_freq_power,'NumChannels',2,...
    'SamplesPerFrame',NumSamples);
% noise with 1/f^aplha spectral characteristic over its entire freq range.
% Ex: alpha = 1 (pink noise), alpha = 2 (brownian noise).

x = step(cn) ; % x: NumSamples x NumChannels
% impose the standard deviation
x = x.*repmat(std_noise./std(x), NumSamples, 1) ; 
% if multiply whole signal by beta, std(x) becomes beta times the old one 
% (var: beta^2 times the old one)

signals_electrodes = struct('Chan_1', x(:,1),'Chan_2', x(:,2)) ; 
signals_electrodes = struct('Chan_1_2', x', 'Chan_2_1', (x(:,[2,1]))') ; 
x_vec = [0:NumSamples-1]./f_s ; 
n_colors = NumChannels ; 
 
if do_plots
    fig_sig = figure('units','normalized','outerposition',...
        [0.1 0.3 0.5 0.6]) ; %, 'Name', ['Subj ',num2str(num_subject),' - ',str_ref, mod_name]) ;
    plot_signals_spectrum(signals_electrodes, x_vec, f_s, n_colors) ;
end

% ======================================================================= %
% ----------> Generating a random periodic signal 
% ======================================================================= %
% used to model the time courses of the different factors.
% Periodic signal: only peaks at k*f_SS, if f_SS is = 1/(period of signal)

% https://nl.mathworks.com/help/ident/ref/idinput.html
% Generate a single-channel periodic random Gaussian input signal with a 
% period of f_s/f_SS samples and NumPeriod periods in the signal. 

% ==> First generate the signal using the entire frequency range 
NumChannels = 2;
Period = f_s/f_SS ; % T_SS/T_s  %50; % number of samples in a period
NumPeriod = 6 ;
NumSamples = Period*NumPeriod ; 

pi_sig = gen_random_periodic_sig(Period, NumChannels,...
    NumPeriod, 'randn') ; % one channel per row
% If the "System Identification Toolbox" is installed on the current
% computer, it is used to generate the random periodic signal
%pi_sig = idinput([Period,NumChannels,NumPeriod],'rgs'); pi_sig = pi_sig' ;


% rgs = random gaussian signal
% one channel per col (time samples along rows)
pi_sig = pi_sig.*std_pi_sig./repmat(std(pi_sig,[],2), 1, NumSamples) ; % impose standard deviation

% ==> then specify a passband
num_harmonics = 3 ; 
band_hz = [0, min(f_s/2, f_SS*(num_harmonics+0.5))] ; % up to num_harmonics harmonics 
band_norm = band_hz./(f_s/2) ;  
% normalized so that 1 corresponds to the Nyquist frequency = f_s/2
%pi_sig_band = idinput([Period,NumChannels,NumPeriod],'rgs',band_norm);
%pi_sig_band = pi_sig_band' ;
pi_sig_band = gen_random_periodic_sig(Period, NumChannels,...
    NumPeriod, 'randn', band_norm) ; % one channel per row

pi_sig_band = pi_sig_band.*std_pi_sig./repmat(std(pi_sig_band,[],2), 1,NumSamples) ; 
% impose standard deviation
 

% Same but with 2 times more harmonics
band_hz = [0, min(f_s/2, f_SS*num_harmonics*2)] ;% up to num_harmonics harmonics 
band_norm = band_hz./(f_s/2) ;  
%pi_sig_band2 = idinput([Period,NumChannels,NumPeriod],'rgs',band_norm);
%pi_sig_band2 = pi_sig_band2' ; 
pi_sig_band2 = gen_random_periodic_sig(Period, NumChannels,...
    NumPeriod, 'randn', band_norm) ; % one channel per row
pi_sig_band2 = pi_sig_band2.*std_pi_sig./repmat(std(pi_sig_band2,[],2), 1, NumSamples) ; 


% ==> Or consider a sum-of-sinusoids signal
% Less well-suited because in this case all harmonics have the same
% amplitude (cfr spectrum)
NumSinusoids = 6 ; 
timeVec = [0:NumSamples-1]./f_s ; 
freqSine = [1:NumSinusoids]'.*f_SS ;
sum_sine_sig = gen_sum_sine_random_phase(freqSine, timeVec) ; 

signals_electrodes = struct(['P_rgs_',num2str(num_harmonics*2)], pi_sig_band2, ...%'P_rgs', pi_sig,...
    ['P_rgs_',num2str(num_harmonics)], pi_sig_band, ...
    'SoSines', sum_sine_sig) ; 
x_vec = timeVec ;  
n_colors = NumChannels ; 
 
if do_plots
    fig_pi_sig = figure('units','normalized','outerposition',...
        [0.1 0.3 0.5 0.6], 'Name', 'Random periodic signals') ; %, 'Name', ['Subj ',num2str(num_subject),' - ',str_ref, mod_name]) ;
    
    norm_signals = 0 ; 
    plot_signals_spectrum(signals_electrodes, x_vec, f_s, n_colors,...
        norm_signals) ;
    xlim([0,f_SS*NumSinusoids*2])
end

% ======================================================================= %
% ----------> Considering a mixture of FAM (factor analysis models)
% ======================================================================= %
% .-------------------------.
% | X_k^i = A*Z_k^i + n_k^i |
% .-------------------------.
% i: indicates the class in {1,2}: nociceptive vs non-nociceptive
% k: indicates the time sample
% where A = [A1, A2]: A1 contains the multi-modal spatial patterns
%                     A2 contains the other spatial patterns
% Z_k^i = [Z_k^{a.i}; Z_k^{b.i}]: 
% Z_k^{a.1} = Z_{k-tau}^{a.1} if k>=tau
%           = 0 if k<tau (NO structured activity for class 1 before tau
%                           sec).
% n_k^i: pink noise, with std defined to reach a given SNR

% simu_data = generate_mix_FAM() ; 
% simu_data: structure with fields samples_cls1, samples_cls2,   
%                   mixingMat and location OR fileNameLoc.
%%%%%%%%%%%%%%%%%%%%%%%%%
NFu = 1 ;               %
NFs = 2 ;               %
stdu = 1 ; % multi-modal activity: can be seen as interference for our purpose!
% should also consider variations in stda btw the 2 classes, cfr AM: the
% multi-modal activity might have different amplitudes in both conditions
% and this should not affect the results...
stds_or_N = 2 ;         %
SNR = 20 ;               % in dB
nChan = 19 ;            %
fs = 500 ;              %
fZ = 0.5 ;              %
% nb of samples per period: round(fs/fZ)
nPeriod = 10 ;          %
nHarmonics = 3 ;        %
% 300ms time lag between 2 classes considered
tau = 0 ; %round(0.3*fs) ;   % time lag between the multi-modal activity, in nb 
                        % of samples
                        % tau should be integer!
inv_freq_power = 1 ;    % 1: pink noise, 2: brownian, 0: white
                        % 2 maybe more realistic than 1! cfr with SNR = -10
                        % or -20 (similar to real EEG SS)
                        % 1: not realistic in time domain!! (too erratic...)
                        % Flicker noise is a type of electronic noise with 
                        % a 1/f, or "pink", power spectral density
                        % 
opt_std = 1 ;           % 1 -->  stds_or_N = stds  
                        % 2 --> stds_or_N = stdN
rng_opt = 'Shuffle' ; %4 ; %'shuffle' ;   %  'Default' % or any number                  
%%%%%%%%%%%%%%%%%%%%%%%%%
simu_data = generate_mix_FAM(NFu, NFs, stdu, stds_or_N, SNR, nChan, ...
    fs, fZ, nPeriod, nHarmonics, tau, inv_freq_power, opt_std, rng_opt) ;
% generate_mix_FAM(NFu, NFs, stdu, stds_or_N, SNR, nChan, ...
%    fs, fZ, nPeriod, nHarmonics, tau, inv_freq_power, opt_std, rng_opt, ...
%    N_RI, std_RI)

A = simu_data.mixingMat ; 
x_cls1 = simu_data.samples_cls1 ;
x_cls2 = simu_data.samples_cls2 ;
mixingMat = simu_data.mixingMat ; 


filter_tmp = pinv(mixingMat) ; 
%tmp = size(mixingMat)
source1 = filter_tmp(1,:)*x_cls1 ; 
% source1 is NOT the same as first_source (the true periodic source)
% because the associated filter should cancel the noise sources!!
% (interferences + random noise)
% --> to recover the same spatial pattern A1, the spatial filter w1 can be
% very different than the one computed by inverting the original A1...

%source1 = mixingMat\x_cls1 ; 
%source1 = source1(1,:) ; 

first_source = simu_data.first_source ; 
all_sources1 = simu_data.all_sources1 ; % all specific sources
all_sources2 = simu_data.all_sources2 ; 
% struct array with fields: samples and possibly location OR fileNameLoc.
% one cell in the array for each class
datasets = struct('samples', {x_cls1; x_cls2}) ; 
% create struct array:
% S(1) = struct('samples', randn(2,3)) ; S(2) = struct('samples', randn(1,3)) ; 
% S(2).newF = 'Tt' % creates the field 'newF'
% OR, directly create a nonscalar struct array:
% S2 = struct('samples', {randn(2,3); randn(1,3); randn(1,2)})

% --> plot the simulated signals and their spectrum
% In freq domain: there is N/2 samples from 0 to fs/2 -> (fs/2)/(N/2) =
% fs/N = fs/(nPeriod*round(fs/fZ))~fZ/nPeriod [Hz] btw 2 samples

x_vec = [0:nPeriod*round(fs/fZ)-1]./fs ;
n1_plot = min(3, nChan) ; 
signals_electrodes = struct(['Cls1_',num2str(n1_plot),'chan'], x_cls1(1:n1_plot,:), ...
    ['Cls2_',num2str(n1_plot),'chan'], x_cls2(1:n1_plot,:), 'Sources1', all_sources1, 'Sources2', all_sources2) ; %, 'S1', source1, 'Orig_S1', first_source) ; 
fig_gen_sig = figure('units','normalized','outerposition',...
    [0.1 0.3 0.5 0.6], 'Name', ['(SNR, stda, stdb_N, opt_std)= ',num2str(SNR), ', ', ...
    num2str(stdu), ', ', num2str(stds_or_N), ', ', num2str(opt_std)]) ;
norm_signals = 0 ; add_spectrum = 1 ; 
xlim_val = [0, 2/fZ] ; % 5 SS periods %NaN ; 
xlim_spec = [0,fZ*nHarmonics*4] ; 
%names_yticks = {'$x_1(t)$', '$x_2(t)$', '$s_1(t)$', 's_2(t)'} ; 
ax_orig = plot_signals_spectrum(signals_electrodes, x_vec, fs, n1_plot,...
    norm_signals, add_spectrum, 1, xlim_val, xlim_spec) ; %, 0,NaN, names_yticks) ;
%set(ax_orig(2), 'YLim', [0, 5000])
% plot_signals_spectrum(signals_electrodes, x_vec, f_s, ...
%     n_colors, norm_signals, add_spectrum, consider_colbar, xlim_val, xlim_spec)

% set(fig_gen_sig,'NextPlot','add'); axes;
% h = title({['SNR = ', num2str(SNR)]; ' '}, 'FontSize', ...
%     20, 'Interpreter', 'Latex') ;
%  set(gca,'Visible','off'); set(h,'Visible','on');

% realistic for SNR<=0...
% because: power of pink noise is spread over all frequencies, while it is
% concentrated on the harmonics for the signals...
% --> maybe: could band pass filter the pink noise before defining its
% variance...

% if there is one of the 2 fields "location" OR "fileNameLoc":
% the scalp plots can be made
do_scalp_plot = 1 ; 
if isfield(simu_data, 'fileNameLoc')
    A = A-repmat(mean(A, 1), nChan, 1) ; 
    fig_A1 = figure() ; hold on ; 
    subplot(121)
    my_topoplot(A(:,1),simu_data.fileNameLoc, 'verbose', 'off', 'style','map') ;
    % style: map, both, contour, fill, blank
    %my_topoplot(A_est(:,1),simu_data.fileNameLoc) ;
    subplot(122)
    my_topoplot(A(:,NFu+1),simu_data.fileNameLoc, 'verbose', 'off', 'style', 'both') ;% both = map + contour lines 
    do_scalp_plot = 1 ;     
    % === compare the first spatial pattern with the fft coefficients at 
    %     freq k*fZ for nH_plot k (k=1, 2,..., nH_plot)
    
    nH_plot = 3 ;
    % nb of harmonics to show on the scalp plots of the spectrum
    [ ~, color_matrix] = get_some_nice_colors(1) ; 
    
    [xtf, freqs] = get_spectrum(x_cls1', fs, 1) ; 
    [xtf_x2, freqs_x2] = get_spectrum(x_cls2', fs, 1) ; 
    fig_fft = figure() ; 
    for idx_harmonic = 1:nH_plot   
        curr_freq = idx_harmonic*fZ ; 
        if idx_harmonic==1
            str_title = ['$|F(X1)|$ at $f_{SS}$'] ; 
            str_title2 = ['$|F(X2)|$ at $f_{SS}$'] ; 
        else
            str_title = [num2str(idx_harmonic),'*$f_{SS}$'] ; 
            str_title2 = [num2str(idx_harmonic),'*$f_{SS}$'] ; 
        end
    
        [~,idx_freq] = min(abs(freqs-curr_freq)) ; 
        main_freq_fft = xtf(idx_freq(1),:) ; 
        % standardization
        main_freq_fft = (main_freq_fft-mean(main_freq_fft))./(std(main_freq_fft)) ; 
        %main_freq_fft = main_freq_fft-repmat(mean(main_freq_fft), 1,
        %nChan) ; 

        [~,idx_freq_x2] = min(abs(freqs_x2-curr_freq)) ; 
        main_freq_fft_x2 = xtf_x2(idx_freq_x2(1),:) ; 
        % standardization
        main_freq_fft_x2 = (main_freq_fft_x2-repmat(mean(main_freq_fft_x2), 1, nChan))./...
            (repmat(std(main_freq_fft_x2), 1, nChan)) ; 
        %main_freq_fft_x2 = main_freq_fft_x2-repmat(mean(main_freq_fft_x2), 1,
        %nChan) ; 
        
        % === X1
        h_X1(idx_harmonic) = subplot(2,nH_plot, idx_harmonic) ; 
        my_topoplot(main_freq_fft,simu_data.fileNameLoc, 'verbose', 'off') ; 
        %set(gca,'units','normalized','position',[0 0 1 1]);
        axis tight
        caxis([-2, 2]) 
        % lims for the colormap % TO ADAPT with min and max values on all plots        
        title(str_title, 'Interpreter', 'LateX', 'FontSize', 14)
        
        % === X2
        h_X2(idx_harmonic) = subplot(2,nH_plot, nH_plot+idx_harmonic) ; 
        my_topoplot(main_freq_fft_x2,simu_data.fileNameLoc, 'verbose', 'off', 'colormap', colormap('jet')) ; 
        % only one colormap by FIGURE!!!
        axis tight
        caxis([-2, 2]) 
        title(str_title2, 'Interpreter', 'LateX', 'FontSize', 14)
        %colormap('winter') ; 
        %axis on
    end
    % tEST
    pos_X1 = get(h_X1, 'Position') ; 
    pos_X2 = get(h_X2, 'Position') ; 
    for idx_h=1:nH_plot
        gap_ypos = pos_X1{idx_h}(2) - (pos_X2{idx_h}(4)+ pos_X2{idx_h}(2)) ; 
        
        %pos_X1{idx_h}(2) = pos_X1{idx_h}(2) - gap_ypos ; 
        %pos_X1{idx_h}(4) = pos_X1{idx_h}(4) + gap_ypos ;
        % cannot change height without changing width! becuase axis equal
        % in my_topoplot
        
        %pos_X2{idx_h}(2) = pos_X2{idx_h}(2) - gap_ypos ; 
        %pos_X2{idx_h}(4) = pos_X2{idx_h}(4) + gap_ypos ;        
        
        set(h_X1(idx_h), 'position', pos_X1{idx_h})
        set(h_X2(idx_h), 'position', pos_X2{idx_h})        
        axis tight
    end
    
    figure(fig_fft) ; 
    
    last_pos = pos_X1{nH_plot} ;  
    last_pos_X2 = pos_X2{nH_plot} ;
    mid_ypos_X2 = last_pos_X2(2)+last_pos_X2(4)*0.5 ; 
    mid_ypos_X1 = last_pos(2)+last_pos(4)*0.5 ; 
    hc = colorbar('peer', h_X2(nH_plot), 'Position', [last_pos(1)+last_pos(3)+0.01,  ...
        0.3, 0.02, 0.5]) ; %mid_ypos_X2, 0.02, mid_ypos_X1-mid_ypos_X2]) ;
    
    %hc.Label.String = '(\mu V)';
    hc.Label.FontSize = 14 ; 
    title(hc, {'a.u.'; ' '}, 'FontSize', 14)
    %hc.Label.Interpreter = 'LateX' ; 
    %%%%%%%%%%%%%%%%%%%
    
    
    do_scalp_plot = 1 ; 
elseif isfield(simu_data, 'location')
    fig_A1 = figure() ;  
    my_topoplot(A(:,1),simu_data.location) ; 
    do_scalp_plot = 1 ;     
else
    do_scalp_plot = 0 ;
end
if do_scalp_plot
    min_axis = floor(2*min(min(A(:,1)), min(A(:,NFu+1))))/2 ; 
    max_axis = ceil(2*max(max(A(:,1)), max(A(:,NFu+1))))/2 ;
    
    figure(fig_A1)
    subplot(121) ; 
    %curr_hdl = gca ; 
    %pos=get(curr_hdl,'position'); set(curr_hdl,'position',[0.8*pos(1) 0.8*pos(2) pos(3) pos(4)])    
    axis tight
    caxis([min_axis, max_axis]) ; 
    title('First spatial activity (A_u)')
    
    subplot(122) ; curr_hdl = gca ; 
    axis tight
    caxis([min_axis, max_axis]) ; 
    title('First specific spatial activity (A_s)')
    % Add the colorbar
    last_pos = get(curr_hdl, 'Position') ;  
    hc = colorbar('peer', curr_hdl, 'Position', [last_pos(1)+last_pos(3)+0.01, last_pos(2)+0.15, ...
        0.02, last_pos(4)-0.2]) ;
    hc.Label.String = '(mu V)';  %hc.Label.FontSize = 14 ; %title(hc, {'\muV'; ' '}, 'FontSize', 14)
end


%r_A = rank(A) %8 = NFa + NFb
%r_sig_cls1 = rank(x_cls1) %19 = nChan
%r_sig_cls2 = rank(x_cls2) %19 = nChan

end

function [res] = gen_sum_sine_random_phase(freqSine, timeVec)
% Generate a sum of sine waves at input frequencies with uniform phases 
% in [0, 2*pi].
% In: - freqSine: a vector of frequencies that have to be considered
%     - timeVec: time vector to evaluate the sine waves.

n_time = length(timeVec) ; 
NumSine = length(freqSine) ; 
r_phi = (2*pi).*rand(NumSine,1) ;% random phases

res = sum(sin(2*pi*freqSine*timeVec +repmat(r_phi, 1, n_time)), 1) ; 
% figure()
% subplot(211)
% plot((sin(repmat(2*pi*freqSine*timeVec, NumSine,1)+repmat(r_phi, 1, n_time)))') ; 
% hold on;
% subplot(212)
% plot(res)

% The fourier transform of a sine wave with a given phase is given by:
% x(t) = sin(w0 t + p1) = (exp(j(w0 t+p1)) - exp(-j(w0 t+p1)))/(2j)--> 
% X(w) = \int_{-\infty}^{\infty} x(t)*exp(-jwt) dt
%      = pi*j*(exp(-j p1)*delta(w+w0) - exp(j p1)*delta(w-w0))

% Therefore, a sum of sine waves with the same frequency and random phases
% is still a sine wave at that frequency, with a phase = the sum of all
% phases of its constituants (taking fourier transform of each wave, doing
% the sum and then iFT).

end

