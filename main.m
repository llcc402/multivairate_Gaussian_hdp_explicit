clear
clc

G = 5; % 5 groups
K = 5; % 5 clusters
D = 2; % 2 dimensional data points
N = 500; % 500 observations in each group

[data, mu, mixing, Z] = data_generate(N, K, G, D);

mat_data = permute(data, [1, 3, 2]);
mat_data = reshape(mat_data, [], D, 1);
plot(mat_data(:,1), mat_data(:,2), 'o')
title('The scatter plot of the whole data set')

figure(2)
plot(data(:,1,1), data(:,2,1), 'o')
title('The scatter plot of the first group')

gamma = 2;
alpha = 1;
actN = 100; % assume there are at most 100 clusters
maxIter = 1000;

tic;
[mu_post, Z_post, mixing_post] = hdp(data, gamma, alpha, actN, maxIter);
toc

figure(3)
hold on 
ix = unique(Z_post(1,:));
xlim([min(data(:,1,1)) - .5, max(data(:,1,1)) + .5])
ylim([min(data(:,2,1)) - .5, max(data(:,2,1)) + .5])
for i = 1:length(ix)
    plot(data(Z_post(1,:) == ix(i), 1, 1), data(Z_post(1,:) == ix(i), 2, 1), 'o')
end
title('The posterior of the indicators for group 1')
hold off