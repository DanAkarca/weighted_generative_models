%% demo weight optimised generative model
% example change
clear; clc;
addpath('/Users/da04/Desktop/Post-doc/Research/weighted_gm');
%% run
nnode = 100;
A = zeros(nnode);
D = abs(randn(nnode)); D = D + D';
m = 100;
modeltype = 'matching';
modelvar = {'powerlaw','powerlaw'};
params = [-2 0.5];
epsilon = 1e-6;
weightupdate = 1;
% ... add the constraint into the optimization rather than fixing after (will speed up)
% ... when adding a weight it should be a weak one
[b c] = weight_optimised_generative_model(A,D,m,modeltype,modelvar,params,epsilon,weightupdate);
%% plot
% plot the communicability over time
for i = 1:m;
    A = squeeze(c(:,:,i));
    w(i) = sum(A,'all');
    [~,q(i)] = modularity_und(A>0);
end
% plot the weights over time
figure; 
plot(w);
figure;
plot(q);
% plot the final graph
Afinal = squeeze(c(:,:,end));
figure; imagesc(Afinal); colorbar;
g = graph(Afinal,'upper'); plot(g);
% plot the strength distribution
figure;
histogram(strengths_und(Afinal))
