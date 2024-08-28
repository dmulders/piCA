function [error_bars_fct_SNR,ah1] = plot_error_bars(...
    xvec, dataset, errors, ylim_val, xlabel_val, ...
    ylabel_val, labels, title_val,filename)
% Plot the input data by linking each pt and adding error bars on top.
%
% In:
%   xvec:   must be the same size as dataset, giving the x position of
%           each data point.
%   dataset:    (m,n) matrix of mean values to plot with bars. Each column  
%                will be ploted as a separated line.
%   errors:     (m,n) matrix with the error values to represent with error
%                bars
%   ylim:   ylim to use
%   labels:     labels of each line (n values should be given)
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be
 
% Default parameters
if nargin<9
    filename = 'unknown_img' ;
end
save_fig = 1 ; 

marker_spec = {'*' , 's', 'p', 'x', '^', 'd', 'o', '+'} ; % , '.'{'.'} ;  
% to plot signgificance below the plots
line_specs = {'-*', '-s', '-p', '-x', '-^', '-d', '-o', '-+'} ; % {'-'}
line_opt_error = '-' ; % line type for the error bar
% to plot the curves
taille = 15 ; 
taille_ticks = 18 ;  
MSize = 8 ; %20 ; % 20 when marker = '.', otherwise: 8
add_legend = 1 ; 

n_subplots = 6 ; % defines fraction of plot to indicate significativity
height_h2 = 0.035 ; % height (rel to fig) to plot the data_signif
gap_ah1_ah2 = 0.18 ; % space necessary for xlabel and xticklabels

[m,n] = size(dataset) ; 
color_matrix = jet(5) ; 
max_colors = size(color_matrix,1) ;
cmap = color_matrix(1:min(max_colors,10),:) ;%jet(n_colors) ;

error_bars_fct_SNR = figure('units','normalized',...
    'outerposition',[0.1 0.1 0.45 0.6]);%, 'Name', cond_name) ;

n_col = size(cmap,1) ; 

ah1 = '' ; %        
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the error bars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mark = length(marker_spec) ; 

idx_mark = 1 ;
for idx_data=1:n
    errorbar(xvec(:,idx_data), dataset(:,idx_data), errors(:, idx_data), ...
        line_opt_error, 'MarkerSize', MSize, 'LineWidth', 2) ; %k. % marker_spec{idx_mark}
    hold on ;
    % LineWidth: for the bars! (and for the data point/lines)
    idx_mark = mod(idx_mark, n_mark) + 1 ;
end
hold on ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add the curves with markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = gca ; ax.ColorOrderIndex = 1 ;
idx_line = 1 ; 
n_lines = length(line_specs) ;
for idx_data=1:n
    h(idx_data) = plot(xvec(:,idx_data), dataset(:,idx_data), line_specs{idx_line}, ...
        'LineWidth', 2, 'MarkerSize', MSize, 'MarkerFaceColor', 'none') ; % h(idx_data) = 
    idx_line = mod(idx_line, n_lines) + 1 ; 
end    

if ~isnan(ylim_val)
    ylim(ylim_val)
end
%set(gca, 'YScale', 'log')
set(gca,'YGrid', 'on', 'XGrid','off', 'FontSize', taille_ticks)
xlabel(xlabel_val, 'FontSize', taille, 'Interpreter', 'Latex') ;
ylabel(ylabel_val, 'FontSize', taille, 'Interpreter', 'Latex') ;
%set(gca,'XScale','log') % 'YScale'
if n>1 && add_legend
   lgd = legend(h, labels, 'Location', 'best', 'Orientation','horizontal') ;%'SouthWest') 
   %lgd = legend('Orientation','horizontal','Interpreter','LateX','fontsize',taille_ticks) ; 
    set(lgd,'color','none');
end
title(title_val, 'FontSize', taille, 'Interpreter', 'Latex') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save figure if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_fig
    eps_fig = 0 ; 
    if eps_fig
        saveas(error_bars_fct_SNR, [filename, '.eps'], 'epsc') ;
    else
        saveas(error_bars_fct_SNR, [filename, '.png'], 'png') ;
    end
end

end


