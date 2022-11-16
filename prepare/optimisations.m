%% testing optimisations
% written by danyal akarca
clear; clc;
%% set paths
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data');
load('consensus');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
%% set data to view
type = 1; % type of dti data
threshold = 4; % thresholded consensus
% example data
a = squeeze(consensus.dk.connectivity(type,threshold,:,:));
b = squeeze(consensus.dk.seed(type,:,:));
c = double(b>0);
% set data to plot
data = a;
%% visualise
% set scaling factor 
sf = 2e-3; % 2:sf, 2e-3:stream
h = figure; h.Position = [100 100 650 650];
subplot(2,2,1);
imagesc(data);
box off; xlabel('Node'); ylabel('Node');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
subplot(2,2,2);
g = graph(data);
plot(g,...
    'linewidth',sf*g.Edges.Weight,...
    'markersize',1e-6+sf*strengths_und(data));
axis off;
subplot(2,2,3);
histogram(degrees_und(data),'edgecolor','w');
ylabel('Frequency'); xlabel('Degree, k');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;
subplot(2,2,4);
histogram(strengths_und(data),'edgecolor','w');
ylabel('Frequency'); xlabel('Strength, s');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;
%% update communicability
% set n
n = 0.05;
% set number of reps either way
rep = 15;
% form a string
repvec = [-rep*n:n:0:n:rep*n];
% nreps
nreps = length(repvec);
% nnode
nnode = size(data,1);
% get postive edges
[a b] = find(triu(data));
lind = find(triu(data));
% compute the initial comunicability
S = diag(sum(data,2));
adj = (S^-.5)*data*(S^-.5);
nc = expm(data);
wc = nc(lind);
% get number of edges
nedge = length(a);
% initialise
norm_comm = nan(nedge,nreps,nnode,nnode);
sum_comm = nan(nedge,nreps);
% loop over edges
for edge = 1:nedge;
    % get index of edge
    ind = [a(edge) b(edge)];
    % get its value
    v = data(ind(1),ind(2));
    % get a vector of delta values
    dv = n*v;
    % set the total changes
    reps = [v-rep*dv:dv:v:dv:v+rep*dv];
    % loop over each
    for r = 1:nreps;
        % get original network
        wsynth = data;
        % update it
        wsynth(ind(1),ind(2)) = reps(r);
        wsynth(ind(2),ind(1)) = reps(r);
        % calculate the normalized communicability
        S = diag(sum(wsynth,2));
        adj = (S^-.5)*wsynth*(S^-.5);
        nc = expm(adj);
        norm_comm(edge,r,:,:) = nc;
        % and the summed communicability
        sum_comm(edge,r) = sum(nc,'all');
    end
end
% plot the sum versus the relative
h = figure; h.Position = [100 100 500 450];
% add middle lines
mid = round(nreps/2);
xline(mid,'linewidth',3,'alpha',.3);
yline(mean(sum_comm(:,mid)),'linewidth',3,'alpha',.3);
hold on;
% keep the coef
coeff = [];
% add in change of weight
for edge = 1:nedge;
    % get number of reps
    nreps = size(sum_comm,2);
    % plot scatter
    x = 1:nreps;
    y = sum_comm(edge,:);
    scatter(x,y,1250,y,'.');
    % fit a linear model
    curve = fit(x',y','poly1');
    % keep the first coef
    coeff(edge) = curve.p1;
    hold on;
    colormap(viridis);
end
% set ticks
tick = [1:5:nreps];
% limits
xlim([tick(1) tick(end)]);
% ticks
xticks(tick); xticklabels(100*repvec(tick));
b = gca; b.TickDir = 'out'; b.FontSize = 18; b.FontName = 'Arial';
xlabel('{\delta\itw_{i,j}} (%)'); ylabel('Total communicabiltiy (\Sigma\itC)');
% correlate the initial communicability with the change
% initial communicability
h = figure; h.Position = [100 100 500 450];
scatter(coeff',wc,2000,'.','markeredgecolor',[.5 .5 .5]);
xlabel('\beta'); ylabel('Initial communicabiltiy \itc_{i,j}');
b = gca; b.TickDir = 'out'; b.FontSize = 20; box off;
[r p] = corr(coeff',wc);