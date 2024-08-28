function [X] = standardize_signals(X, only_center)
% Standardize the signals in X, containing one channel per row and one 
% time step per column. 
% There might be a third dimension (eg for epochs).
% We standardize by removing the mean along time and then dividing by 
% standard deviation.

if nargin<2
    only_center = 0 ; 
end

%[n_C, n_T, n_ep] = size(X) ; 
n_T = size(X, 2) ; 

mean_along_T = mean(X, 2) ; 

if only_center
    X = X - repmat(mean_along_T,1,n_T) ; 
else
    std_along_T = std(X, 0, 2) ; % 0: norm by N-1; 1: norm by N
    %[i_C, i_T,i_ep] = find(abs(std_along_T)<10^(-6))
    std_along_T(abs(std_along_T)<10^(-6)) = 1  ;% avoid division by 0 if one 
                                            % std is close to 0.
    X = (X - repmat(mean_along_T,1,n_T))./repmat(std_along_T,1,n_T) ; 
end

end

