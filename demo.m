%% demo weight optimised generative model
% example change
clear; clc;
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
%% load data
% load the data
load('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat');
%% run weighted generative model
% target
type = 1;
threshold = 6;
Wtgt = squeeze(consensus.dk.connectivity(type,threshold,:,:));
Wseed = squeeze(consensus.dk.seed(type,:,:)); Wseed(Wseed<0.5)=0;
D = consensus.dk.euclidean;
% compute observed statistics
x = cell(4,1);
x{1} = sum(weight_conversion(Wtgt,'normalize'),2);
x{2} = clustering_coef_wu(weight_conversion(Wtgt,'normalize'));
x{3} = betweenness_wei(weight_conversion(Wtgt,'normalize'));
x{4} = D(triu(Wtgt,1) > 0);
% binary generative model parameters
nnode = 68;
seed = zeros(nnode,nnode);
m = nnz(Wtgt)/2;
modeltype = 'matching';
modelvar = {'powerlaw','powerlaw'};
params = [-2 0.35];
epsilon = 1e-6;
% weighted genreative model parameters
weighted_model = struct;
weighted_model.update = 1; % flag for weighted genreative models optimisation
weighted_model.omega = 1; % the parameterisation of communicability when optimised
weighted_model.start = 1; % when to start optimisation
weighted_model.optimisation.resolution = 0.05; % how many weight alternations are sampled before the gradient is taken
weighted_model.optimisation.samples = 2; % how many samples are taken to infer the gradient
weighted_model.optimisation.alpha = 0.02; % the update coefficient
% run the weighted genreative model
[bb cc dd] = weight_optimised_generative_model(seed,D,m,modeltype,modelvar,params,epsilon,...
    weighted_model);
% get the final network
wfinal = dd(:,:,end);
% options string
opt = sprintf('omega=%g, start=%g, resolution=%g, samples=%g, alpha=%g',...
    weighted_model.omega,...
    weighted_model.start,...
    weighted_model.optimisation.resolution,...
    weighted_model.optimisation.samples,...
    weighted_model.optimisation.alpha);
%% do an evaluation
% iniialise
K = [];
y = cell(4,1);
% compute on the simulation
y{1} = sum(weight_conversion(wfinal,'normalize'),2);
y{2} = clustering_coef_wu(weight_conversion(wfinal,'normalize'));
y{3} = betweenness_wei(weight_conversion(wfinal,'normalize'));
y{4} = D(triu(wfinal,1) > 0);
for j = 1:4
    K(j) = fcn_ks(x{j},y{j});
end
E = max(K,[],2);
% labels
ks_labels = string({'s','wc','wb','d'});
% visualise
h = figure; h.Position = [10 10 1500 300];
for i = 1:4;
    subplot(1,5,i); 
    histogram(x{i},10,'edgecolor','w'); hold on; histogram(y{i},10,'edgecolor','w');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
    xlabel(ks_labels(i)); ylabel('Frequency');
end
hold on;
subplot(1,5,5);
bar(K,'edgecolor','w'); xticklabels(ks_labels); ylim([0 1]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
ylabel('KS'); 
% title
sgtitle(opt);
%% plot outcome
% come some statistical measures over time
nstat = 2;
statistics = zeros(m,nstat);
for step = 1:m;
    % take network
    w = squeeze(dd(:,:,step));
    % compute statistics
    statistics(step,1) = sum(w,'all'); % total weight
    [~,statistics(step,2)] = modularity_und(w); % modularity
end
% plot the weights over time
h = figure; h.Position = [10 10 800 300];
subplot(1,2,1);
plot(statistics(:,1));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time'); ylabel('Network weight'); box off;
subplot(1,2,2);
plot(statistics(:,2));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time');  ylabel('Modularity Q'); box off;
% title
sgtitle(opt);

% plot the final graph
h = figure; h.Position = [10 10 800 300];
subplot(1,2,1);
imagesc(weight_conversion(Wtgt,'normalize'));
c = colorbar; c.Label.String = 'Weight'; colormap(viridis);
xlabel('Node'); ylabel('Node');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;
subplot(1,2,2);
imagesc(weight_conversion(wfinal,'normalize'));
c = colorbar; c.Label.String = 'Weight'; colormap(viridis);
xlabel('Node'); ylabel('Node');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;
% title
sgtitle(opt);
%% ks function
function kstat = fcn_ks(x1,x2)
    binEdges    =  [-inf ; sort([x1;x2]) ; inf];
    binCounts1  =  histc (x1 , binEdges, 1);
    binCounts2  =  histc (x2 , binEdges, 1);
    sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
    sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);
    sampleCDF1  =  sumCounts1(1:end-1);
    sampleCDF2  =  sumCounts2(1:end-1);
    deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
    kstat = max(deltaCDF);
end