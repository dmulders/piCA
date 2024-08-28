function [W,D] = sort_W_based_on_D(W, D, order)
% Sort the columns of W and associated diagonal entries in D in decending
% order of the real value of the diagonal elements in D (D is a matrix,
% not vector).
%
% Inputs:
%   - W: matrix with as many columns as D.
%   - D: diagonal matrix (! not vector) according to which the columns of W
%           will be sorted.
%   - order: 'descend' (default) or 'ascend' indicates the criterion for
%               the sort.
% ----------------------------------------------------------------------- %

% Dounia Mulders - dounia.mulders@uclouvain.be

if nargin<3
    order = 'descend' ; 
    % order: 'ascend', 'descend'
end

D_vec = diag(D) ;
[~, idx_sorted] = sort(real(D_vec), order) ; 
D = diag(D_vec(idx_sorted)) ; 
W = W(:,idx_sorted) ; 

end

