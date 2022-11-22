%% demo weight optimised generative model
% example change
clear; clc;
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/models/');
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
Wseed = squeeze(consensus.dk.seed(type,threshold,:,:));
D = consensus.dk.euclidean;
% compute observed binary statistics
Atgt = double(Wtgt>0);
xb = cell(4,1);
xb{1} = degrees_und(Atgt)';
xb{2} = clustering_coef_bu(Atgt);
xb{3} = betweenness_bin(Atgt)';
xb{4} = D(triu(Atgt,1) > 0);
% compute observed weighted statistics
xw = cell(4,1);
xw{1} = sum(weight_conversion(Wtgt,'normalize'),2);
xw{2} = clustering_coef_wu(weight_conversion(Wtgt,'normalize'));
xw{3} = betweenness_wei(weight_conversion(Wtgt,'normalize'));
xw{4} = D(triu(Wtgt,1) > 0);
% binary generative model parameters
nnode = 68;
% seed = Wseed;
% seed = double(Wseed>0);
seed = zeros(nnode);
m = nnz(Wtgt)/2;
mseed = nnz(seed)/2;
modeltype = 'matching';
modelvar = {'powerlaw','powerlaw'};
params = [-2.5 0.35];
epsilon = 1e-6;
% weighted genreative model parameters
weighted_model = struct;
weighted_model.update = 1; % flag for weighted genreative models optimisation
weighted_model.start = 1; % when to start optimisation (this is on top of the initial seed network connectivity)
weighted_model.optimisation.function = 1; % communicability, 1. normalised communicability,2
weighted_model.optimisation.resolution = 0.05; % how many weight alternations are sampled before the gradient is taken
weighted_model.optimisation.samples = 2; % how many samples are taken to infer the gradient
weighted_model.optimisation.alpha = .02; % the update coefficient [must be a smaller range]
% run the weighted genreative model
[bb cc dd] = weight_optimised_generative_model(seed,D,m,modeltype,modelvar,params,epsilon,...
    weighted_model);
% get the final network
wfinal = dd(:,:,end);
% options string
opt = sprintf('start=%g, resolution=%g, samples=%g, alpha=%g,func=%g',...
    weighted_model.start,...
    weighted_model.optimisation.resolution,...
    weighted_model.optimisation.samples,...
    weighted_model.optimisation.alpha,...
    weighted_model.optimisation.function);
%% do an evaluation
% compute compute binary
Abin = double(wfinal>0);
Kb = [];
yb = cell(4,1);
yb{1} = degrees_und(Abin)';
yb{2} = clustering_coef_bu(Abin);
yb{3} = betweenness_bin(Abin)';
yb{4} = D(triu(Abin,1) > 0); % can cause an error if there are negative weights!
for j = 1:4
    Kb(j) = fcn_ks(xb{j},yb{j});
end
Eb = max(Kb,[],2);
% compute weighted
Kw = [];
yw = cell(4,1);
yw{1} = sum(weight_conversion(wfinal,'normalize'),2);
yw{2} = clustering_coef_wu(weight_conversion(wfinal,'normalize'));
yw{3} = betweenness_wei(weight_conversion(wfinal,'normalize'));
yw{4} = D(triu(wfinal,1) > 0); % can cause an error if there are negative weights!
for j = 1:4
    Kw(j) = fcn_ks(xw{j},yw{j});
end
Ew = max(Kw,[],2);
% labels
ks_labels = string({'s','wc','wb','d'});
% visualise
h = figure; h.Position = [10 10 1500 300];
for i = 1:4;
    subplot(1,5,i); 
    histogram(xw{i},10,'edgecolor','w'); hold on; histogram(yw{i},10,'edgecolor','w');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
    xlabel(ks_labels(i)); ylabel('Frequency');
end
hold on;
subplot(1,5,5);
bar(Kw,'edgecolor','w'); xticklabels(ks_labels); ylim([0 1]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
ylabel('KS'); 
% title
sgtitle(opt);
%% plot outcome
% come some statistical measures over time
nstat = 3;
statistics = zeros(m-mseed,nstat);;
for step = 1:m-mseed;
    % take network
    w = squeeze(dd(:,:,step));
    % compute statistics
    statistics(step,1) = sum(w,'all'); % total weight
    [~,statistics(step,2)] = modularity_und(w); % modularity
    statistics(step,3) = efficiency_wei(w);
end
% plot the weights over time
h = figure; h.Position = [10 10 800 300];
subplot(1,3,1);
plot(statistics(:,1));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time'); ylabel('Network weight'); box off;
subplot(1,3,2);
plot(statistics(:,2));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time');  ylabel('Modularity Q'); box off;
subplot(1,3,3);
plot(statistics(:,3));
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