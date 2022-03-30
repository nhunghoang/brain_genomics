function U1 = diffusion_embedding(W, c, r)
% This function computes a diffusion embedding of an input timeseries or
% correlation matrix. It largely follows the Marguiles et al. (2016) approach.
% inputs:
%       W,      timeseries or correlation matrix
%       c,      number of output gradients
%       r,      number of computed singular vectors (to speed up)
%
% output:
%       U1,     gradient maps


% data dimensions
[n, t] = size(W);

% number of output gradient
if ~exist('c', 'var')
    c = 3;
end

% number of singular vectors
if ~exist('r', 'var')
    r = 2 * c;
end

% process timeseries
if n ~= t                           % if timeseries
    X = norm_normalize(W, 2);       % normalize to unit norm
    W = X * X';                     % compute cosine similarity
end

% sparsify
W(W < prctile(W, 90, 2)) = 0;       % keep 10% strongest out-connections for each node
W = norm_normalize(W, 2);           % normalize to unit norm

% compute sparse svd
[U0, S0, ~] = svds(W, r);

% correlation of correlations
W = U0 * (S0^2) * U0';              % compute consine similarity of cosine similarities
W = max(W, 0);                      % remove negative correlations (there should be few)
K = W; clear W                      % store in variable K and free up memory

% diffusion probability
sqrt_k = sqrt(sum(K, 2));
K = K ./ (sqrt_k * sqrt_k');        % compute diffusion kernel
P = K; clear K;                     % store in variable P and free up memory


% gradients
[U1, ~, ~] = svds(P, c + 1);        % spatial singular vectors

% postprocess
U1 = U1 ./ sqrt_k;                  % convert to eigenvectors
U1 = normalize(U1(:, 2:end), 1);    % normalize and remove first (constant) eigenvector

end


function W = norm_normalize(W, dim)
% normalization by L2 norm
W = W ./ sqrt(sum(W.^2, dim, 'omitnan'));
W(isnan(W)) = 0;

end
