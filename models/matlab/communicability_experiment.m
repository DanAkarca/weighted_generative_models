% generate multiple networks and examine communicability optimisations
% this phenomena seems to be related to how modular the network is
% written by danyal akarca
clear; clc;
% add paths
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
% load calm data
load('/imaging/astle/users/da04/PhD/qsiprep_data/dti/calm_qsiprep.mat');
%% take a subject's data
% set subject
sub = 6;
A_calm = squeeze(calm_qsiprep.aal116.connectivity(sub,1,:,:));
% set nnode
nnode = 116;
% visualise
g = graph(A_calm);
figure; 
plot(g,...
    'LineWidth',.0005*g.Edges.Weight,...
    'MarkerSize',1e-3+.0005*sum(A_calm,2));
axis off;
% plot the degree distribution
figure;
subplot(1,2,1);
histogram(degrees_und(A_calm));
xlabel('Degree, k');
ylabel('Frequency');
% plot the strength distribution
subplot(1,2,2);
histogram(strengths_und(A_calm));
xlabel('Strength, s');
ylabel('Frequency');
%% take the average calm
% average over calm
A_calm_av = squeeze(mean(squeeze(calm_qsiprep.brainnetome246.connectivity(:,1,:,:)),1));
% visualise
g = graph(A_calm_av);
figure; 
plot(g,...
    'LineWidth',.0005*g.Edges.Weight,...
    'MarkerSize',1e-3+.0005*sum(A_calm_av,2));
axis off;
% plot the degree distribution
figure;
subplot(1,2,1);
histogram(degrees_und(A_calm_av));
xlabel('Degree, k');
ylabel('Frequency');
% plot the strength distribution
subplot(1,2,2);
histogram(strengths_und(A_calm_av));
xlabel('Strength, s');
ylabel('Frequency');
%% compute random density matched networks
% do the MS algorithm 
Arand = randmio_und(A_calm_av,10);
% visualise
g = graph(Arand);
figure; 
plot(g,...
    'LineWidth',.0005*g.Edges.Weight,...
    'MarkerSize',1e-3+.0005*sum(Arand,2));
axis off;
% plot the degree distribution
figure;
subplot(1,2,1);
histogram(degrees_und(Arand));
xlabel('Degree, k');
ylabel('Frequency');
% plot the strength distribution
subplot(1,2,2);
histogram(strengths_und(Arand));
xlabel('Strength, s');
ylabel('Frequency');
%% compute total communicability as a function of weight changes
% set the data to plot
dataplot = A_calm;
% set the extent of the weight changes
dx = [0:0.5:2];
ntweak = length(dx);
% get the indicies of edges to change
u = triu(ones(nnode),1);
% find their location
[x y] = find(u);
nedge = length(x);
% initialise
output = zeros(nedge,ntweak);
% loop through edges
for edge = 1:nedge;
    % loop through iterative changes to each edge
    for tweak = 1:length(dx);
        % get a copy
        b = dataplot;
        % update the edge value
        b(x(edge),y(edge)) = dx(tweak)*dataplot(x(edge),y(edge));
        % compute the total communicability
        s = diag(sum(b,2));
        u = s^-.5 *b *s^-.5;
        output(edge,tweak) = sum(expm(u),'all');
    end
end
% plot
h = figure; h.Position = [100 100 600 400];
for edge = 1:nedge;
    plot(output(edge,:));
    hold on;
end
b = gca; b.TickDir = 'out'; b.FontSize = 18; b.FontName = 'Arial';
xlabel('{\delta\itw_{i,j}} (%)'); ylabel('Total communicabiltiy (\Sigma\itC)');
box off;
%% evaluate where the difference occurs between communicability and shortest paths
% set data to plot
data = normalize(A_calm,'range');;
% calculate the shortest path length between each node pair
[spl,hops,Pmat] = distance_wei_floyd(data,'inv');
% calculate the communicability
s = diag(sum(data,2));
u = s^-.5 * data *s^-.5;
c = expm(u);
% correlate the two measures
figure;
scatter(-log10(c(:)),log10(spl(:)),'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('-log_{10}(Communicability)'); ylabel('log_{10}(Shortest path)');
% color nodes by the type of connection
%% compute the distortion of signals
% set the data
data = A_calm; % can normalize by strength here
% calculate the shortest path length between each node pair
[spl,hops,Pmat] = distance_wei_floyd(data,'inv');
% compute the communicability 
s = diag(sum(data,2));
u = s^-.5 * data *s^-.5;
c = expm(u);
% distortion is defined as the probability that a walker will not take the shortest path
% follow the shortest path and calculate a cumulative probability it was
% taken by the diffusion operator
[u v] = find(data);
% keep the distortion
distortion = zeros(nnode);
% loop over all node pairs
for i = 1:length(u);
    % get the source
    source = u(i);
    for j = 1:length(v);
        % get the target
        target = v(j);
        % retrieve the shortest path
        pth = retrieve_shortest_path(source,target,hops,Pmat);
        % initialise the specific path
        distortion_path = [];
        % follow the path
        for k = 1:length(pth)-1;
            % take the shortest path steps
            hop_a = pth(k);
            hop_b = pth(k+1);
            % calculate the probability (relative to all others) that the signal diffused through this path
            % take the diffusion from that node
            diffuse = c(hop_a,:);
            % take the distortion at this hop
            short_path = diffuse(hop_b);
            distortion_path(k) = 1-(diffuse(hop_b)/sum(diffuse));
        end
        % multiply all hops (i.e., the probability it did all hops)
        distortion(u(i),v(j)) = prod(distortion_path,'all');
    end
end
%% use this distortion measure
h = figure; h.Position = [100 100 500 400];
imagesc(distortion);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12;
box off; l = colorbar; l.Label.String = 'Distortion';
xlabel('Node'); ylabel('Node');
% correlate with the original measures
h = figure; h.Position = [100 100 1000 300];
subplot(1,3,1);
scatter(distortion(:),c(:),'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Distortion'); ylabel('Communicability');
subplot(1,3,2);
scatter(distortion(:),spl(:),'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Distortion'); ylabel('Shortest path length');
subplot(1,3,3);
scatter(c(:),spl(:),'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Communicability'); ylabel('Shortest path length');