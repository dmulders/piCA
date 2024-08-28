function [ ] = test_piCA(varargin)
% Analyze the performances of the piCA algorithms implemented in
% piCA_compute on simulated EEG data.
% The data are generated using generate_mixed_sig (see test_EEG_model for an 
% example). 
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be

close all ; 

reload_res = false ;        % reload previously saved simulations  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 500 ;                  % Sampling frequency
f_SS = 0.5 ;                % Frequency of activity to highlight
                            % --> nb of samples per period: round(fs/f_SS)
NF = 1 ;                    % Number of hidden factors
stds_or_N = 1 ;             % SD of the hidden source or noise
opt_std = 1 ;               % 1 -->  stds_or_N = stds  
                            % 2 --> stds_or_N = stdN
nChan = 2 ;                %
nPeriod = 10 ;              %
nHarmonics = 5 ;            % nb of harmonics for the factors Z
inv_freq_power = 1 ; %1 ;   % inverse frequency power for the PSD of the 
                            % additive noise.
   
nH = 10 ;                   % number of harmonics to consider, v in {3,4}. 
                            % 2*nH signals will be  
                            % considered for the reference periodic signal.
freqs_vec_to_max = ([-nH:-1, 1:nH]).*f_SS ;
                            % For version==4.
comp_fft_pattern = 1 ;      % compare the results to the FFT spatial 
                            % pattern at f_SS (canonical angle wrto the
                            % first spatial pattern is computed).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
n_runs = 10 ;               % nb of random generations of the signals
all_SNR = -30:2:0 ;         % 31 values
                            % not realistic to take large positive SNR
                            % SNR (dB) defines var of additive noise 
n_SNR = length(all_SNR) ;   %   
plot_degree = 1 ;           % whether to show the angles in degrees or rad
show_spatial_patterns = 1 ; % only if nChan = 2 

n_harm_save = 10 ;          % To compute Mpi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
version_piCA = [1, 2, 3, 4] ; % in {1,2, 3, 4}

color_matrix = jet(5) ; 

max_colors = size(color_matrix,1) ;
cmap = color_matrix(1:min(max_colors,10),:) ;%jet(n_colors) ;
lwdt = 2 ;  

str_std_other = '_stds' ; 
if opt_std==2
    str_std_other = '_stdN' ; 
end
str_all_std = [str_std_other, num2str(stds_or_N)] ; 

% params for all methods
gen_params = [num2str(n_runs), 'runs_', num2str(n_SNR), 'SNR_FP', ...
    num2str(inv_freq_power), str_all_std, '_', num2str(nChan), ...
        'C_F', num2str(NF), '_', num2str(nPeriod),'T'] ;
% FP = frequency power 
% C = channels

n_methods = length(version_piCA) ; 

% *** Data structures for the results
angles_all_methods = NaN*ones(n_runs, n_SNR,n_methods) ; 
% angle between the first col of A and A_est
first_angles_all_methods = NaN*ones(n_runs, n_SNR, n_methods) ; 
periodicity_methods = NaN*ones(n_runs, n_SNR, n_methods) ; 
xcorr_S1 = NaN*ones(n_runs, n_SNR, n_methods) ; 

if comp_fft_pattern
    first_angles_fft_topo = NaN*ones(n_runs, n_SNR) ; 
    all_fft_topo = zeros(nChan, n_runs, n_SNR) ; 
end

str_caption_methods = cell(n_methods, 1) ; 

for idx_method = 1:n_methods
    piCA_v = version_piCA(idx_method) ; 
    disp(['===== Considering method ', num2str(piCA_v), '(',...
        num2str(idx_method),'/',num2str(n_methods),')'])

    if ~reload_res
        all_A = cell(n_runs, n_SNR) ; 
        all_W = cell(n_runs, n_SNR) ; 
        all_A_est = cell(n_runs, n_SNR) ; 
        % angle between the subspaces spanned by A and A_est
        all_angles = NaN*ones(n_runs, n_SNR) ; 
        % angle between the first col of A and A_est
        all_first_angles = NaN*ones(n_runs, n_SNR) ; 

        all_Mpi = NaN*ones(n_runs, n_SNR) ; 
        all_xcorr = NaN*ones(n_runs, n_SNR) ; 

        for idx_SNR = 1:n_SNR
            curr_SNR = all_SNR(idx_SNR) ; 
            disp(['Considering SNR = ', num2str(curr_SNR), '(',...
                num2str(idx_SNR),'/',num2str(n_SNR),')'])
            for idx_run=1:n_runs
                % index in an array (n_SNR,n_runs)
                idx_run_SNR = (idx_SNR-1)*n_runs+idx_run ; 
                rng_opt = idx_run_SNR ; 
                % --> each method is applied to the same signals
                % === Generate the data
                simu_data = generate_mixed_sig(NF, stds_or_N, curr_SNR, nChan, ...
                    fs, f_SS, nPeriod, nHarmonics, inv_freq_power, opt_std, rng_opt) ; 

                A = simu_data.mixingMat ; 
                samples = standardize_signals(simu_data.samples,0) ;
                all_A{idx_run, idx_SNR} = A ; 
                first_source = simu_data.first_source ;
                
                % === Compute the spatial patterns using LPDA
                T_SS_nb = round(fs/f_SS) ; 
                datasets = struct('samples', samples, 'period', T_SS_nb) ;   
                [W, D] = piCA_compute(datasets, nH, piCA_v, ...
                    freqs_vec_to_max, fs) ; 

                A_est = pinv(W') ;
                if ~isreal(A_est)
                    warning(['Estimated A is not real, only keep real part'])
                    A_est = real(A_est) ; 
                end

                n_filters = size(W,2) ; 
                % normalize the spatial patterns
                A_est = A_est./(repmat(sqrt(sum(A_est.^2, 1)), nChan,1)) ; 

                all_A_est{idx_run, idx_SNR} = A_est ; 
                all_W{idx_run, idx_SNR} = W ;   

                best_filter = W(:,1) ; 
                if ~isreal(W)
                    W = real(W) ;
                    if ~isreal(best_filter)
                        disp(['-- optimal filter isn''t real: '])                              
                        best_filter = W(:,1) ; 
                    end  
                end

                NF_used = NF ; 
                if NF>1
                    NF_used = floor(NF/2) ; 
                end
                NF_used = min(NF_used, n_filters) ;

                all_angles(idx_run, idx_SNR) = subspace(A(:,1:NF_used), ...
                    A_est(:, 1:NF_used)) ; 
                all_first_angles(idx_run, idx_SNR) = subspace(A(:,1), ...
                    A_est(:, 1)) ;
                best_filtered_sig = best_filter'*samples ; 
                % ATTENTION, for method 2, x_cls1 can be averaged over periods!!!
                all_Mpi(idx_run, idx_SNR) = compute_Mpi(best_filtered_sig, ...
                    fs, f_SS, n_harm_save) ;
                if length(best_filtered_sig)~=length(first_source)
                    warning(['! first_source has length ', num2str(length(first_source)), ...
                        ' whereas best_filtered_sig has length ', ...
                        num2str(length(best_filtered_sig))])
                else
                    all_xcorr(idx_run, idx_SNR) = abs(corr(best_filtered_sig',...
                        reshape(first_source, length(first_source), 1))) ;
                end

                if comp_fft_pattern && idx_method==1
                    % Compute and store FFT topo
                    topo_tmp = compute_FFT_topo(samples', fs, f_SS) ; 
                    topo_tmp = topo_tmp./(sqrt(sum(topo_tmp.^2))) ; 
                    first_angles_fft_topo(idx_run, idx_SNR) = ...
                        subspace(A(:,1), topo_tmp) ; 

                    all_fft_topo(:,idx_run, idx_SNR) = topo_tmp ; 
                end
            end
        end
    end

    % === Filenames ===
    str_method = ['meth', num2str(piCA_v)] ;
    if any(piCA_v==[3,4])
        str_method = [str_method, '_nH', ...
            num2str(nH)] ; 
    end
    filename_dir = ['./results/', str_method, '/'] ;
    filename = [filename_dir, 'data_', gen_params] ; 

    dir_fft = ['./results/FFT_topo/'] ;
    filename_fft = [dir_fft, 'data_',gen_params] ; 

    if reload_res
        % --> loading the data
        all_curr_data = load([filename, '.mat']) ; 
        all_angles = all_curr_data.all_angles ; 
        all_first_angles = all_curr_data.all_first_angles ; 
        all_Mpi = all_curr_data.all_Mpi ; 
        all_xcorr = all_curr_data.all_xcorr ; 
        all_A_est = all_curr_data.all_A_est ; 
        all_A = all_curr_data.all_A ; 
        all_W = all_curr_data.all_W ; 

        if comp_fft_pattern && idx_method==1
            % Load FFT topo                
            all_fft_data = load([filename_fft, '.mat']) ;
            first_angles_fft_topo = all_fft_data.first_angles_fft_topo ; 
            all_fft_topo = all_fft_data.all_fft_topo ; 
        end
    else
        % --> saving the data
        if ~exist(filename_dir, 'dir') % 0 or 7 if it exists
            mkdir(filename_dir)
        end
        %save(filename) % save the whole workspace
        save([filename, '.mat'], 'all_A', 'all_A_est', 'NF', 'stds_or_N', ...
            'nChan', 'fs', 'f_SS', 'nPeriod', 'nHarmonics', ...
            'inv_freq_power', 'opt_std', 'all_angles', 'all_first_angles', ...
            'all_W', 'all_Mpi', 'all_xcorr') ;   

        if comp_fft_pattern && idx_method==1
            if ~exist(dir_fft, 'dir') % 0 or 7 if it exists
                mkdir(dir_fft)
            end
            save([filename_fft, '.mat'], 'first_angles_fft_topo', 'all_fft_topo') ; 
        end
    end

    str_caption_methods{idx_method} = ['\piCA', ...
        ' v', num2str(piCA_v)] ; 
    angles_all_methods(:,:, idx_method) = all_angles ; 
    first_angles_all_methods(:,:, idx_method) = all_first_angles ; 
    periodicity_methods(:,:,idx_method) = all_Mpi ;  
    xcorr_S1(:,:,idx_method) = all_xcorr ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check: display spatial patterns   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------------------------------- %   
    if show_spatial_patterns
        % show the true and estimated spatial patterns (only if nChan in
        % [2, 3])
        if nChan==2
            figure('Name', str_caption_methods{idx_method}, 'units','normalized',...
                'outerposition',[0.1 0.1 0.45 0.6]) ; % ['Method ',num2str(method)]
            hold on ; 
            set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

            show_W_plot = 1 ; % show NF_used filters

            A = all_A{n_runs, n_SNR} ;
            A_est = all_A_est{n_runs, n_SNR} ;   
            W = all_W{n_runs, n_SNR} ;   
            n_filters = min(size(A,2), size(A_est,2)) ; 
            % normalize the spatial patterns
            A_est = A_est./(repmat(sqrt(sum(A_est.^2, 1)), nChan,1)) ; 

            NF_used = min(NF, n_filters) ;              
            n_remaining_F = n_filters-NF_used ; 
            % Consider vectors in 2 first quadrants only (through central
            % symetry if necessary)
            A = give_equil_col_2first_quadrants(A) ;
            A_est = give_equil_col_2first_quadrants(A_est) ;
            W = give_equil_col_2first_quadrants(W) ;

            % -- true patterns
            h_lgd(1:NF_used) = plot([zeros(1, NF_used); A(1,1:NF_used)], ...
                [zeros(1, NF_used); A(2,1:NF_used)], ...
                '-', 'LineWidth', lwdt) ; hold on ;   
            if show_W_plot
                lgd_str = cell(1, NF_used*2 + n_remaining_F + NF_used) ; 
            else
                lgd_str = cell(1, NF_used*2 + n_remaining_F) ; 
            end
            lgd_str(1:NF_used) = {'True Fs'} ; 

            % -- estimated patterns (only Fs)
            ax = gca ; ax.ColorOrderIndex = 1 ;
            h_est_lgd(1:NF_used) = plot([zeros(1, NF_used); A_est(1,1:NF_used)], ...
                [zeros(1,NF_used); A_est(2,1:NF_used)], ...
                '-.', 'LineWidth', lwdt) ; hold on ;   
            lgd_str((NF_used+1):(2*NF_used)) = {'Est. Fs'} ; 
            % Attention: nb of cols of A_est can be smaller
            % than nb of cols A (=NF)
            if n_remaining_F>0
                h_est_lgd((NF_used+1):(NF_used+n_remaining_F)) = plot([zeros(1, n_remaining_F); ...
                    A_est(1,(NF_used+1):(NF_used+n_remaining_F))], ...
                    [zeros(1,n_remaining_F); A_est(2,(NF_used+1):(NF_used+n_remaining_F))], ...
                    '-.', 'LineWidth', lwdt) ; hold on ; 
                lgd_str((2*NF_used+1):(2*NF_used+n_remaining_F)) = {'Remaining. F'} ; 
            end

            if show_W_plot
                h_est_lgd((NF_used+n_remaining_F+1):(2*NF_used+n_remaining_F)) = plot([zeros(1, NF_used); W(1,1:NF_used)],...
                    [zeros(1, NF_used); W(2,1:NF_used)], ...
                    '--', 'LineWidth', lwdt) ; hold on ; 
                lgd_str((2*NF_used+n_remaining_F+1):...
                    (2*NF_used+n_remaining_F+NF_used)) = {'W used'} ; 
            end

            if comp_fft_pattern
                % show it with ALL methods, or can aslso impose a
                % method number here
                fft_tmp = all_fft_topo(:, n_runs, n_SNR) ; 
                h_est_lgd(2*NF_used+n_remaining_F+1) = plot([0; fft_tmp(1)], ...
                [0; fft_tmp(2)], 'k-.', 'LineWidth', lwdt) ; hold on ; 
                lgd_str(end+1) = {'FFT'} ; % 'FFT topo'
            end

            grid on ; 
            legend([h_lgd, h_est_lgd], lgd_str) ; 
        elseif nChan==3
            % Plot 3D vectors if nChan=3
        end
    end
    % ------------------------------------------------------------------- %    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_folder = ['./results/'] ; 
if ~exist(str_folder, 'dir') % 0 or 7 if it exists
    mkdir(str_folder)
end
taille = 20 ; 
% --- canonical angles between the true and found subspaces
str_angle = 'angle' ; 
str_params = [num2str(n_methods), 'meth_', gen_params,...
    '_nH', num2str(nH)] ;
filename_all_angles = [str_folder,str_angle,'_', str_params] ;
%line_specs = {'-o', '-s', '-p', '-x', '-^', '-d'} ; % {'-'}
line_specs = {'-*', '-s', '-p', '-x', '-^', '-d', '-o', '-+'} ; % {'-'}

if plot_degree
    angles_all_methods = angles_all_methods.*(180/pi) ; 
    first_angles_all_methods = first_angles_all_methods.*(180/pi) ; 
    str_angle = '($^{\circ}$)' ; 
    ylim_val = [0, 90] ; 
    if comp_fft_pattern
        first_angles_fft_topo = first_angles_fft_topo.*(180/pi) ; 
    end
else
    str_angle = '(rad)' ; 
    ylim_val = [0, pi/2] ; 
end

xvec = repmat(all_SNR', 1, n_methods) ; 
% ensure that the datasets have the right dimensions, even if there is only
% one method --> reshape
dataset = reshape(squeeze(mean(angles_all_methods, 1)),n_SNR,n_methods) ; 
% each column of datasets: a different method
errors = reshape(squeeze(std(angles_all_methods, 1, 1)),n_SNR,n_methods) ;

% Angle between all the true and estimated patterns
plot_error_bars(xvec, dataset, errors, ylim_val, 'SNR (dB)', ...
    ['Canonical angle ', str_angle], str_caption_methods, ...
    '', filename_all_angles) ;


% --- canonical angles between the FIRST true and found pattern
str_first_angle = 'angle1' ; 
filename_first_angles = [str_folder,str_first_angle,'_', str_params] ;

xvec_angle1 = xvec ; str_caption_angle1 = str_caption_methods ; 
dataset = reshape(squeeze(mean(first_angles_all_methods, 1)),n_SNR,n_methods) ; 
% each column of datasets: a different method
errors = reshape(squeeze(std(first_angles_all_methods, 1, 1)),n_SNR,n_methods) ;
if comp_fft_pattern
    xvec_angle1 = repmat(all_SNR', 1, n_methods+1) ; 
    dataset = [dataset, (mean(first_angles_fft_topo, 1))'] ; 
    errors = [errors, (std(first_angles_fft_topo, 1, 1))'] ; 
    str_caption_angle1(end+1) = {'FFT'} ; % 'FFT topo'
    % first_angles_fft_topo: (n_runs, n_SNR)
    % first_angles_all_methods: NaN*ones(n_runs, n_SNR,n_methods) ; 
    
    %%%%%%%%%%%
    n_sets = length(str_caption_angle1) ; 
    color_matrix = jet(10) ; 
    max_colors = size(color_matrix,1) ;
    cmap = color_matrix(1:min(max_colors,n_sets),:) ;%jet(n_colors) ;
    idx_FFT = mod(n_sets-1, length(cmap)) + 1 ; 
    cmap(idx_FFT,:) = [0,0,0] ; 
    %%%%%%%%%%%%%%%%%%%
end

plot_error_bars(xvec_angle1, dataset, errors, ylim_val, 'SNR (dB)', ...
    ['First canonical angle ', str_angle], str_caption_angle1, ...
    '', filename_first_angles) ;

% --- Periodicity Mpi of the first extracted sources
str_Mpi = 'Mpi' ; 
filename_Mpi = [str_folder, str_Mpi,'_', str_params] ;
% Mpi for the first estimated source
dataset = reshape(squeeze(mean(periodicity_methods, 1)),n_SNR,n_methods) ; % each column of datasets: a different method
errors = reshape(squeeze(std(periodicity_methods, 1, 1)),n_SNR,n_methods) ; 
tmp_data = dataset - errors ; 
ylim_val(1) = floor(min(tmp_data(:))) ; 
tmp_data = dataset + errors ;
ylim_val(2) = ceil(max(tmp_data(:))) ;

plot_error_bars(xvec, dataset, errors, ylim_val, 'SNR (dB)', ...
    ['$M_{\pi}(s_1^1)$'], str_caption_methods, ...
    '', filename_Mpi) ;  

% --- Correlation between true and estimated first source Zs1
str_xcorr = 'xcorr' ; 
filename_xcorr = [str_folder, str_xcorr,'_', str_params] ;
% Mpi for the first estimated source
dataset = reshape(squeeze(mean(xcorr_S1, 1)),n_SNR,n_methods) ; % each column of dataset: a different method
errors = reshape(squeeze(std(xcorr_S1, 1, 1)),n_SNR,n_methods) ; 
tmp_data = dataset - errors ; 
ylim_val(1) = 0 ; %floor(min(tmp_data(:))*20)/20 ; 
tmp_data = dataset + errors ;
ylim_val(2) = 1 ; %ceil(max(tmp_data(:))*20)/20 ;

plot_error_bars(xvec, dataset, errors, ylim_val, 'SNR (dB)', ...
    ['Corr($s_1^1(t), \mathbf{z}_{s,1}^{1}(t)$)'], str_caption_methods, ...
    '', filename_xcorr) ; 

end

function  [res_mat] = give_equil_col_2first_quadrants(in_mat)
% Return a matrix of the same dimensions as in_mat, where each column
% vector is adapted so that it lies in the 2 first quadrants ([0, pi]).
% in_mat should have a first dimension of 2 or 3 (otherwise we will not
% display its columns as vectors anyway).

[n_rows,n_col] = size(in_mat) ; 
res_mat = in_mat ; 
if any(n_rows==[2, 3]) 
    for idx_col = 1:n_col
        curr_col = in_mat(:,idx_col) ; 
        if curr_col(2)<0
            % when negative y component, do a central symmetry
            % NOT an axial symmetry wrto the x-axis, because with the
            % central symmetry we stay on the same direction in space.
            res_mat(:, idx_col) = curr_col.*(-1) ; 
        end
    end
end

end



