% Demo weighted generative model code
% Written by Danyal Akarca, 2023

% clear workspace and command window
clear; clc;
% **** set directory **** set to the weighted generative model code directory
dir = '/Users/da04/Desktop/code_share/';

% go to directory
cd(dir);
% add code to path
addpath(genpath(dir));
% load
load demo_data;

% get empirical data
Wtgt = demo_data.Wtgt;
seed = demo_data.seed;
D = demo_data.D;
coordinates = demo_data.coordinates;

% set model parameters
modeltype = 'matching';
modelvar = {'powerlaw','powerlaw'}; 

%% demo 1a: run a single exemplar weighted generative model

% set how many connections to generate
m = 80; % e.g., 80 connections can generate a small weighted network and does it within ~30-45 seconds

% set an example parameter combination
eta = -3.2; gamma = 0.38; alpha = 0.05; omega = 0.9;
params = [eta gamma alpha omega];

% set the settings
weighted_model = struct;
weighted_model.update = 1; % runs the weighted algorithm
weighted_model.start = 1; % starts optimising as connections form
weighted_model.optimisation.function = 2; % normalised communicability
weighted_model.optimisation.resolution = 0.05; % default
weighted_model.optimisation.samples = 5; 

% run the generative model
[bb,cc,dd] = weighted_generative_model(seed,D,m,modeltype,modelvar,params,1e-5,weighted_model);

%% demo 1b: visualise this single run

% get the final network
Wsim = squeeze(dd(:,:,end));

% visualise
h = figure;
h.Position = [100 100 800 750];
g = graph(Wsim);
plot(g,...
    'XData',coordinates(:,1),...
    'YData',coordinates(:,2),...
    'ZData',coordinates(:,3),...
    'MarkerSize',5*strengths_und(Wsim)+1,...
    'NodeColor',[.3 .3 .3],...
    'EdgeColor',[.7 .7 .7],...
    'LineWidth',4*g.Edges.Weight+1e-5,...
    'NodeLabel',[]);
axis off;

%% demo 1c: visualise the evolution of this network

h = figure;
h.Position = [100 100 600 500];
for i = 1:m;
    x = squeeze(dd(:,:,i));
    s = strengths_und(x)+1;
    s = ceil(s);
    g = graph(x);
        z = plot(g,...
        'XData',coordinates(:,1),...
        'YData',coordinates(:,2),...
        'ZData',coordinates(:,3),...
        'MarkerSize',3*s,...
        'NodeColor',[183 201 226]./255,...
        'EdgeColor',[.9 .9 .9],...
        'LineWidth',4*g.Edges.Weight+1e-5,...
        'EdgeAlpha',.8,...
        'NodeLabel',[]);
    axis off;
    view(i+10,40);
    pause(0.4);
    hold on;
    set(gcf,'color','w');
end

%% demo 2: evaluate example weighted generative models

% set wiring model parameters
nparams = 2; % set number of combinations to run, this can be altered
eta = linspace(-3.7,-2.7,nparams);
gamma = linspace(0.35,0.45,nparams);
alpha = linspace(0.02,0.1,nparams);
omega = linspace(0.8,1.01,nparams);
params = [eta' gamma' alpha' omega'];
epsilon = 1e-6;

% weighted genreative model parameters
weighted_model = struct;
weighted_model.update = 1;
weighted_model.start = 1;
weighted_model.optimisation.function = 2;
weighted_model.optimisation.resolution = 0.05;
weighted_model.optimisation.samples = 5; 

% run and evaluate the weighted generative model
[output] = evaluate_weighted_generative_model(Wtgt,seed,D,modeltype,modelvar,params,epsilon,weighted_model);
