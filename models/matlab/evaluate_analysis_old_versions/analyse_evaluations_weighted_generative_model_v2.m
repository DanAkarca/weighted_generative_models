%% analyse weighted model outputs (v2 - note this for variable omega, but not varied parameters)
% written by danyal akarca, university of cambridge, 2022

% clear
clear; clc;

% add paths
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/models/matlab');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');

% load empirical data
load '/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat' consensus;

% change directory
directory = '/imaging/astle/users/da04/Postdoc/weighted_gm/model_outputs_131222b';
cd(directory);

%% load all networks

% get list
list = dir(fullfile(directory,'*.mat'));
list = string({list.name});
nnet = length(list);
% initialise key outputs
alpha = zeros(nnet,1);
run = zeros(nnet,1);
energy = zeros(nnet,2);
nks = 4;
ks = zeros(nnet,nks,2);
ks_labels = string({'Strength','Clustering','Betweenness','Edge length'});
% load all networks
for net = 1:nnet;
    % load the network
    load(list(net));
    % extract model related variables
    alpha(net) = output.model.optimisation.alpha;
	omega(net) = output.model.optimisation.omega;
    run(net) = output.run;
    % extract evaluations
    energy(net,:) = [output.evaluation.binary.energy output.evaluation.weighted.energy];
    ks(net,:,1) = [output.evaluation.binary.ks]; ks(net,:,2) = [output.evaluation.weighted.ks];
    % clear output
    clear output
    % display
    disp(sprintf('Network %g of %g loaded',net,nnet));
end
% get indices for networks over the same omega, alpha and run
[~,vo,io] = unique(omega);
[~,va,ia] = unique(alpha);
[~,vr,ir] = unique(run+1);
% number of unique alpha and runs
nomega = length(vo);
nalpha = length(va);
nrun = length(vr);

%% plot energy over the whole sample
% visualise
h = figure; h.Position = [100 100 1000 400];
subplot(1,2,1);
histogram(energy(:,1),'edgecolor','w'); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
subplot(1,2,2);
histogram(energy(:,2),'edgecolor','w'); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;

%% plot model fits over omega and alpha

% visualise energy relationship to alpha
h = figure; h.Position = [100 100 500 400];
scatter(alpha,energy(:,2),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
xlabel('\alpha'); ylabel('Energy');
ylim([0 1]);

% visualise ks relationship to alpha
h = figure; h.Position = [100 100 1000 800];
for k = 1:nks;
    subplot(2,nks/2,k);
    scatter(alpha,ks(:,k,2),200,'.');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
    xlabel('\alpha'); ylabel('KS'); 
    title(ks_labels(k),'FontWeight','normal','FontSize',25,'FontName','Arial');
    ylim([0 1]);
end

% visualise energy relationship to omega
h = figure; h.Position = [100 100 500 400];
scatter(omega,energy(:,2),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
xlabel('\omega'); ylabel('Energy');
ylim([0 1]);

% visualise ks relationship to omega
h = figure; h.Position = [100 100 1000 800];
for k = 1:nks;
    subplot(2,nks/2,k);
    scatter(omega,ks(:,k,2),200,'.');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25; box off;
    xlabel('\omega'); ylabel('KS'); 
    title(ks_labels(k),'FontWeight','normal','FontSize',25,'FontName','Arial');
    ylim([0 1]);
end

%% plot the energy and ks landscapes

h = figure; h.Position = [100 100 600 500];
scatter(omega,alpha,10000,energy(:,2),'.'); colorbar;
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\omega'); ylabel('\alpha');

%% compare distributions to the empirical

% plot an example empirical comparison
% find the network with the lowest energy
[i v] = min(energy);

% select a network
network = v(2);
% load final weighted simulation
load(list(network));
Wsynth = squeeze(output.network.weighted(:,:,end));
% load the empirical
Wtgt = squeeze(consensus.dk.connectivity(1,6,:,:));
% define the comparison
x1 = strengths_und(weight_conversion(Wtgt,'normalize'))';
x2 = strengths_und(weight_conversion(Wsynth,'normalize'))';
% take common bins
binEdges =  [-inf;sort([x1;x2]);inf];
% compute the bin counts
binCounts1 = histc (x1,binEdges,1);
binCounts2 = histc (x2,binEdges,1);
% compute a cumulative probability
sumCounts1 = cumsum(binCounts1)./sum(binCounts1);
sumCounts2 = cumsum(binCounts2)./sum(binCounts2);
% compute the cumulative density
sampleCDF1 = sumCounts1(1:end-1);
sampleCDF2 = sumCounts2(1:end-1);
% compute the diffference 
deltaCDF = abs(sampleCDF1-sampleCDF2);
% compute the largest difference (KS)
kstat = max(deltaCDF);

% visualise the CDF
h = figure; h.Position = [100 100 500 400];
plot(sampleCDF1,'linewidth',2);
hold on;
plot(sampleCDF2,'linewidth',2);
hold on;
yline(kstat,'--','linewidth',2);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 25;
box off;
xlabel('X'); ylabel('Cumulative probability');

% visualise the plots

% plot the empirical statistics
y = cell(4,1);
y{1} = sum(weight_conversion(Wsynth,'normalize'),2);
y{2} = real(clustering_coef_wu(weight_conversion(Wsynth,'normalize')));
y{3} = betweenness_wei(weight_conversion(Wsynth,'normalize'));
D = consensus.dk.euclidean;
y{4} = D(triu(Wsynth,1) > 0);

% labels
ks_labels = string({'s','wc','wb','d'});
% visualise
h = figure; h.Position = [10 10 1500 300];
for i = 1:4;
    subplot(1,4,i); 
    histogram(y{i},10,'edgecolor','w');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
    xlabel(ks_labels(i)); ylabel('Frequency');
end

%% visualise the relationship between KS
figure; 
imagesc(corr(ks(:,:,2)));

%% extract networks at specific parameters

