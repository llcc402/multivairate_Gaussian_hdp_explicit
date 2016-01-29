% Input:
%     N        a scalar. The number of observations in each group.
%     K        a scalar. The number of clusters in the whole data set.
%     G        a scalar. The number of groups.
%     D        a scalar. The number of dimension.
% Output:
%     data     a matrix of order N * D * G. 
%     mu       a matrix of order K * D. The i-th row is the center of
%              cluster i.
%     mixing   a matrix of order G * K. The i-th row is the mixing measure
%              of group i. 
%     Z        a matrix of order G * N. The element Z(i,j) is the cluster
%              number of point j in cluser i. 
function [data, mu, mixing, Z] = data_generate(N, K, G, D)
if nargin < 1
    N = 500;
end
if nargin < 2
    K = 5;
end
if nargin < 3
    G = 5;
end
if nargin < 4
    D = 2;
end
%--------------------------------------------------------------------------
% STEP 1: Generate mixing measure
%--------------------------------------------------------------------------
mixing = zeros(G, K);
for k = 1:G
    n = randi(K); % number of clusters in group k
    ix = randperm(K, n); % the selected clusters
    mixing(k,ix) = rand(1, n);
    mixing(k,:) = mixing(k,:) / sum(mixing(k,:));
end

%--------------------------------------------------------------------------
% STEP 2: Generate latent variables
%--------------------------------------------------------------------------
Z = zeros(G, N);
for i = 1:G
    [~,~,Z(i,:)] = histcounts(rand(1, N), [0, cumsum(mixing(i,:))]);
end

%--------------------------------------------------------------------------
% STEP 3: Generate cluster centers
%-------------------------------------------------------------------------- 
sigma = eye(D) * (4 * D)^2;
mu = mvnrnd(zeros(1,D), sigma, K); % one row per center.

%--------------------------------------------------------------------------
% STEP 3: Generate the data set
%--------------------------------------------------------------------------
data = zeros(N, D, G);
for i = 1:G
    for j = 1:N
        data(j,:,i) = mvnrnd(mu(Z(i, j),:), eye(D));
    end
end

end