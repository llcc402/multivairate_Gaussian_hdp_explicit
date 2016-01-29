function [mu_post, Z_post, mixing_post] = hdp(data, gamma, alpha, actN, maxIter)
if nargin < 4
    actN = 100;
end
if nargin < 5
    maxIter = 1000;
end

%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------

[N, D, G] = size(data); % N: the number of data points in each group.      
                        % D: the dimension of the data points
                        % G: the number of the groups

Z_post = ones(G, N);

M_data = permute(data, [1, 3, 2]);
M_data = reshape(M_data, [], D, 1); % reshape the data into a 2d mat.

sigma = std(M_data); % the std of the data set
mu_post = mvnrnd(zeros(1,D), eye(D) .* diag(sigma), actN);
mixing_post = zeros(G, actN);

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    
    % sample G0
    a = histcounts(Z_post(:), 1:actN+1);
    b = [cumsum(a(2:end), 'reverse'), 0];
    a = a + 1;
    b = b + gamma;
    V = betarnd(a,b);
    G0 = V;
    V = cumprod(1-V);
    G0(2:end) = G0(2:end) .* V(1:end-1);
    
    % sample G1:GM
    for i = 1:G
        counts = histcounts(Z_post(i,:), 1:actN+1);
        a = alpha * G0 + counts;
        b = [cumsum(a(2:end), 'reverse'), 0];
        V = betarnd(a, b);
        mixing_post(i,:) = V;
        V = cumprod(1 - V);
        mixing_post(i,2:end) = mixing_post(i,2:end) .* V(1:end-1);
    end
    
    % sample mu
    id = reshape(Z_post', [], 1); % reshape the indicators to a vector. The 
                                 % order of the vector is group1, group2,
                                 % ...
    
    ix = accumarray(id, 1:length(id), [], @(x){x}); % ix is a cell
    for i = 1:length(ix)
        if ~isempty(ix{i})
            mu_post(i,:) = mean(M_data(ix{i}, :));
        end
    end
    ix_0 = setdiff(1:actN, unique(Z_post(:)));
    mu_post(ix_0, :) = mvnrnd(zeros(1, D), ...
                              eye(D) .* diag(sigma),...
                              length(ix_0));
    
    % sample Z_post
    for i = 1:size(Z_post, 1)
        for j = 1:size(Z_post, 2)
            log_prob = log(mixing_post(i,:)') - ...
                       sum((mu_post - repmat(data(j,:,i), actN, 1)) .^ 2, 2);
            log_prob = log_prob - max(log_prob);
            prob = exp(log_prob);
            [~, ~, Z_post(i,j)] = histcounts(rand(1), [0; cumsum(prob)]);
        end
    end
    
end

end
