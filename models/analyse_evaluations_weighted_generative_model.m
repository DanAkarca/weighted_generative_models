%% analyse weighted model outputs
% written by danyal akarca, university of cambridge, 2022

clear; clc;
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/models/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
directory = '/imaging/astle/users/da04/Postdoc/weighted_gm/model_outputs_221122';
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
ks = zeros(nnet,4,2);
% load all networks
for net = 1:nnet;
    % load the network
    load(list(net));
    % extract model related variables
    alpha(net) = output.model.optimisation.alpha;
    run(net) = output.run;
    % extract evaluations
    energy(net,:) = [output.evaluation.binary.energy output.evaluation.weighted.energy];
    ks(net,:,1) = [output.evaluation.binary.ks]; ks(net,:,2) = [output.evaluation.weighted.ks];
    % clear output
    clear output
    % display
end
% get indices for networks over the same alpha and run
[~,va,ia] = unique(alpha);
[~,vr,ir] = unique(run+1);
% number of unique alpha and runs
nalpha = length(va);
nrun = length(vr);

% question: why is the binary model giving consistent energy values? check that Abin is calculated correctly (e.g., negative weights)

%% plot energy over the whole sample


%% plot model fits at each alpha

% take a summary energy for all alpha values evaluated
energy_alpha = zeros(nalpha,2);
for alpha = 1:nalpha;
    % take the mean energy for each alpha
    energy_alpha(alpha,:) = mean(energy(ia(va==alpha),:),1);
end

%% plot model fits over runs within each alpha


%% extract networks at specific parameters

